#!/usr/bin/perl
$sim_name = "geant";
$refit = 1;
$fake_err = 0;
$npar = 7;
@param_label = qw(b_t0 b_a0 b_a1 b_a2 b_b0 b_b1 b_b2);
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
$param_file = "input_params_theta_${sim_name}_$fit_detector.dat";
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
	printf "  %15s =  %10.4g +/- %10.4g\n", $param_label[$p], $par[$p], $dpar[$p];
	}
printf "  Chisq = %.2f Number of Fit Parameters = $ndf\n",$chisq;
printf "  Reduced chisq = %.3f\n", $redchisq;

# Write a file with the fitted parameters
$param_outfile = "parameters_theta_${sim_name}_${fit_detector}_${data_set}.dat";

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
	my ($val, $x, $p);
	my ($B1, $B2, $B3, $B4, $B5, $B6, $B7);
	$x = $_[0]; $p = $_[1]; 

	$B1 = $par[0] - $par[1]*$par[4];
	$B2 = $par[1];
	$B3 = $par[2];
	$B4 = $par[3];
	$B5 = -($par[2]*$par[4] + $par[1]*$par[5]);
	$B6 = -($par[3]*$par[4] + $par[2]*$par[5] + $par[1]*$par[6]);
	$B7 = -($par[3]*$par[5] + $par[2]*$par[6]);
	$B8 = -$par[3]*$par[6];

	$val = $B1 +$B2*$p +$B3*$p*$x +$B4*$p*$x**2
		+$B5*$x +$B6*$x**2 +$B7*$x**3 +$B8*$x**4;
	return $val;
	}
sub function_err_2
	{
	my ($val, $x, $p, $dx, $dp, $ddx, $ddp);
	$x = $_[0]; $p = $_[1]; 
	$dx = $_[2]; $dp = $_[3]; 

	#$B1 = $par[0] - $par[1]*$par[4];
	$B2 = $par[1];
	$B3 = $par[2];
	$B4 = $par[3];
	$B5 = -($par[2]*$par[4] + $par[1]*$par[5]);
	$B6 = -($par[3]*$par[4] + $par[2]*$par[5] + $par[1]*$par[6]);
	$B7 = -($par[3]*$par[5] + $par[2]*$par[6]);
	$B8 = -$par[3]*$par[6];

	$ddx = $B3*$p + 2.*$B4*$p*$x + $B5 +2.*$B6*$x + 3.*$B7*$x**2 + 4.*$B8*$x**3;
	$ddp = $B2 +$B3*$x +$B4*$x**2;
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
		$term = $th_in[$i] - function($x_out[$i],$th_out[$i]);
		$term = $term**2;
		$term = $term / function_err_2($x_out[$i],$th_out[$i],$d_x_out[$i],$d_th_out[$i]);
		$sum += $term;
		}
	return $sum;
	}
