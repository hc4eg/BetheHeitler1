#!/usr/bin/perl
$sim_name = "geant";
$refit = 1;
$fake_err = 0;
$npar = 10;
@param_label = qw(a_0 a_x a_2x a_3x a_th a_2th a_3th a_xth a_2xth a_x2th);
$central_energy = 30.;

if($#ARGV >= 0) {$fit_detector = $ARGV[0];}
else {$fit_detector = 1;}
if($#ARGV >= 1) {$data_set = $ARGV[1];}
else {$data_set = "ideal";}
print "Fitting the energy for detector $fit_detector\n";

$input_file = "${sim_name}_rays_${data_set}_${fit_detector}.dat";
print "Input data file name: $input_file\n";

# First read in the data
open( IN, $input_file) || die "Cannot open file $input_file\n";
$n = -1;
while(<IN>)
	{
	next if /^#/;
	chomp;
	@data = split(/\s+/);
	if($data[0] =~ /Event:/)
		{ # a new event
		if(!$det_OK && $n != -1)
			{ print "WARNING: No data from detector $fit_detector in last event.\n";}
		else    {
			$n++;
			$ev[$n] = $data[1];
			}
		$det_OK = 0;
		}
	elsif ($data[0] =~ /Input:/)
		{
		$energy_in[$n] = $data[1];
		$de_in[$n] = $data[2];
		$x_in[$n] = $data[3];
		$y_in[$n] = $data[4];
		$th_in[$n] = $data[5];
		$ph_in[$n] = $data[6];
		}
	elsif ($data[0] =~ /Detector:/)
		{
		$detector = $data[1];
		if($detector == $fit_detector)
			{ $det_OK = 1;} else {$det_OK = 0;}
		}
	elsif ($data[0] =~ /VDC:/ && $det_OK)
		{
		$x_out[$n] = $data[1];
		$d_x_out[$n] = 0.3;
		$y_out[$n] = $data[2];
		$d_y_out[$n] = 0.9;
		$th_out[$n] = $data[3];
		$d_th_out[$n] = 18.;
		$ph_out[$n] = $data[4];
		$d_ph_out[$n] = 13.;
		}
	}
close(IN);
$num = $n + 1;
print "Input rays = $num\n";

#for($i = 0; $i < $num; $i++)
#	{
#	print "Event $ev[$i]:\n";
#	print "energy_in[$i] = $energy_in[$i] ";
#	print "delta_in[$i] = $de_in[$i] ";
#	print "x_in[$i] = $x_in[$i] ";
#	print "y_in[$i] = $y_in[$i] ";
#	print "theta_in[$i] = $th_in[$i] ";
#	print "phi_in[$i] = $ph_in[$i] ";
#	print "\n";
#	print "x_out[$i] = $x_out[$i] ";
#	print "y_out[$i] = $y_out[$i] ";
#	print "theta_out[$i] = $th_out[$i] ";
#	print "phi_out[$i] = $ph_out[$i] ";
#	print "\n";
#	}

# read a file with initial guesses
$param_file = "input_params_energy_${sim_name}_$fit_detector.dat";
open( IN, $param_file) || die "Cannot open file $param_file\n";
for($i = 0; $i < $npar; $i++)
	{ $par[$i] = 0; $step[$i] = 0; }

while(<IN>)
	{
	chomp;
	@data = split(/\s+/);
	for($i = 0; $i < $npar; $i++)
		{
		if($data[0] eq $param_label[$i])
			{
			$par[$i] = $data[1];
			$step[$i] = $data[2];
			$pmin[$i] = $data[3];
			$pmax[$i] = $data[4];
			}
		}
	}
close(IN);

print "Input parameters.... from \"$param_file\"\n";
for($i = 0; $i < $npar; $i++)
	{
	printf "  %10s =  %10.4g step %10.4g min %10.3g max %10.3g\n", $param_label[$i], $par[$i], $step[$i], $pmin[$i], $pmax[$i];
	}
# do the fit
$chisq = search();
print "First fit parameters....\n";
for($i = 0; $i < $npar; $i++)
	{
	printf "  %10s =  %10.4g step %10.4g dpar %10.4g\n", $param_label[$i], $par[$i], $step[$i], $dpar[$i];
	}
printf "Chisq = %.4g\n",$chisq;

# refit with better steps
$ndf = 0;
for($i = 0; $i < $npar; $i++)
	{
	if($step[$i] > 0)
		{
		$ndf++;
		if($fake_err) { $step[$i] /= 5.;}
		else	{
			if($dpar[$i] < $step[$i]/20.)
				{ $step[$i] = $step[$i]/20.;}
			else
				{ $step[$i] = $dpar[$i]/2.;}
			}
		}
	}
if($refit) { $chisq = search(); }

if($num > $ndf)
	{ $redchisq = $chisq/($num-$ndf);}
else
	{ $redchisq = 0;}
	
print "Fitted parameters....\n";
for($p = 0; $p < $npar; $p++)
	{
	printf "  $param_label[$p]  =  %.4g +/- %.4g\n", $par[$p], $dpar[$p];
	}
printf "  Chisq = %.2f Number of Fit Parameters = $ndf\n",$chisq;
printf "  Reduced chisq = %.3f\n", $redchisq;

# Write a file with the fitted parameters
$param_outfile = "parameters_energy_${sim_name}_${fit_detector}_${data_set}.dat";
open(OUT, ">$param_outfile") || die "Cannot open file $outfile\n";
for($p = 0; $p < $npar; $p++)
	{
	printf OUT "$param_label[$p] %.4g %.4g\n", $par[$p], $dpar[$p];
	}
close(OUT);
print "Parameters written to: $param_outfile\n";


# Write a file with the fitted function
#$outfile = "fit.dat";
#open(OUT, ">$outfile") || die "Cannot open file $outfile\n";
#for($i = 0; $i < $num; $i++)
#	{
#	$fit = function($x_out[$i], $th_out[$i] );
#	print OUT "$de_in[$i] $fit\n";
#	}
#close(OUT);

exit;

sub function
	{
	my ($val, $x, $p);
	$x = $_[0]; $p = $_[1]; 
	$val = $par[0]
		+ $par[1]*$x  + $par[2]*$x*$x + $par[3]*$x**3
		+ $par[4]*$p + $par[5]*$p*$p + $par[6]*$p**3
		+ $par[7]*$x*$p + $par[8]*$x*$x*$p + $par[9]*$x*$p*$p;
	return $val;
	}
sub function_err_2
	{
	my ($val, $x, $p, $dx, $dp, $ddx, $ddp);
	$x = $_[0]; $p = $_[1]; 
	$dx = $_[2]; $dp = $_[3]; 
	$ddx = $par[1] + 2.*$par[2]*$x + 3.*$par[3]*$x**2
		+ $par[7]*$p + $par[8]*2.*$x*$p +$par[9]*$p*$p;
	$ddp = $par[4] + 2.*$par[5]*$p + 3.*$par[6]*$p**2
		+ $par[7]*$x + $par[8]*$x*$x +$par[9]*2.*$p*$x;
	$val = ($ddx*$dx)**2 + ($ddp*$dp)**2;
	#print "function_err_2 x = $x p = $p dx = $dx dp = $dp val = $val\n";
	return $val;
	}
	


sub search
	{
	my($better, $chi, $chiadd, $chisub, $i, $p);
	$chi = chisq();
	my $num_steps = 0;
	do
	    {
	    $better = 0;
	    for($p=0; $p<$npar; $p++)
	    	{
	    	if($step[$p] != 0.)
	    		{
	    		#print "par[$p] = $par[$p], step[$p] = $step[$p]\n";
	    		$par[$p] += $step[$p];
			if($par[$p] > $pmax[$p]) { $chiadd = $chi; }
			else { $chiadd = chisq(); }
	    		$par[$p] -= 2*$step[$p];
	    		if($par[$p] < $pmin[$p]) {$chisub = $chi; }
			else { $chisub = chisq(); }
	    		#print "chi = $chi, chiadd = $chiadd, chisub = $chisub\n";
	    		if($chisub < $chi)
	    			{
	    			$chi = $chisub;
	    			$better = 1;
	    			}
	    		elsif ($chiadd < $chi)
	    			{
	    			$chi = $chiadd;
	    			$better = 1;
	    			$par[$p] += 2*$step[$p];
	    			}
	    		else
	    			{
	    			$par[$p] += $step[$p];
	    			}
	    		}
	    	}
	    $num_steps++;
	    } while($better);
	print "Converged after $num_steps steps\n";
	# Calculate the error in the parameters
	for($p=0; $p<$npar; $p++)
		{
		if($step[$p] != 0)
			{
			$par[$p] += $step[$p];
			$chiadd = chisq();
			$par[$p] -= 2*$step[$p];
			$chisub = chisq();
			$par[$p] += $step[$p];
			$dpar[$p] = $step[$p] / sqrt( ($chiadd + $chisub)/2. - $chi);
			}
		else
			{
			$dpar[$p] = 0.;
			}
		}
	return $chi;
	}

sub chisq
	{
	my ($sum, $i, $term);
	$sum = 0;
	for($i=0;$i<$num;$i++)
		{
		$term = $de_in[$i] - function($x_out[$i],$th_out[$i]);
		$term = $term**2;
		$term = $term / function_err_2($x_out[$i],$th_out[$i],$d_x_out[$i],$d_th_out[$i]);
		$sum += $term;
		}
	return $sum;
	}
