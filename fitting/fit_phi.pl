#!/usr/bin/perl
$sim_name = "geant";
$refit = 1;
$fake_err = 0;
$npar = 6;
@param_label = qw(c_0 c_1 c_2 c_3 c_4 c_5);

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
# read a file with initial guesses
$param_file = "input_params_phi_${sim_name}_$fit_detector.dat";
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
			last;
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
			if($dpar[$i] < 20.*$step[$i])
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
	printf "  %15s =  %10.4g +/- %10.4g\n", $param_label[$p], $par[$p], $dpar[$p];
	}
printf "  Chisq = %.2f Number of Fit Parameters = $ndf\n",$chisq;
printf "  Reduced chisq = %.3f\n", $redchisq;

# Write a file with the fitted parameters
$param_outfile = "parameters_phi_${sim_name}_${fit_detector}_${data_set}.dat";

open(OUT, ">$param_outfile") || die "Cannot open file $outfile\n";
for($p = 0; $p < $npar; $p++)
	{
	printf OUT "$param_label[$p] %.7g %.7g\n", $par[$p], $dpar[$p];
	}
close(OUT);
print "Parameters written to: $param_outfile\n";

# Write a file with the fitted function
#$outfile = "fit.dat";
#open(OUT, ">$outfile") || die "Cannot open file $outfile\n";
#for($i = 0; $i < $num; $i++)
#	{
#	$fit = function($x_out[$i], $th_out[$i] );
#	print OUT "$th_in[$i] $fit\n";
#	}
#close(OUT);

exit;

sub function
	{
	my ($val, $x, $y, $t, $p);
	$x = $_[0]; $y = $_[1];
	$t = $_[2]; $p = $_[3];

	$val = $par[0]
		+ $par[1]*$p
		+ $par[2]*$p**2
		+ $par[3]*$t
		+ $par[4]*$x
		+ $par[5]*$y;
	return $val;
	}
sub function_err_2
	{
	my ($val, $x, $y, $t, $p, $dx, $dy, $dt, $dp);
	my ($ddx, $ddy, $ddt, $ddp);
	$x = $_[0]; $y = $_[1];
	$t = $_[2]; $p = $_[3]; 
	$dx = $_[4]; $dy = $_[5];
	$dt = $_[6]; $dp = $_[7]; 

	$ddx = $par[4];
	$ddy = $par[5];
	$ddt = $par[3];
	$ddp = $par[1] + 2.*$par[2];
	$val = ($ddx*$dx)**2
		+ ($ddy*$dy)**2
		+ ($ddt*$dt)**2
		+ ($ddp*$dp)**2;
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
		$term = $ph_in[$i] - function($x_out[$i],$y_out[$i], $th_out[$i], $ph_out[$i]);
		$term = $term**2;
		$term = $term / function_err_2($x_out[$i],$y_out[$i], $th_out[$i], $ph_out[$i],
				$d_x_out[$i],$d_y_out[$i], $d_th_out[$i], $d_ph_out[$i]);
		$sum += $term;
		}
	return $sum;
	}
