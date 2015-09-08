#!/usr/bin/perl
use Math::Trig;

$sim_name = "geant";
if($#ARGV >= 0) {$fit_detector = $ARGV[0];}
else {$fit_detector = 1;}
if($#ARGV >= 1) {$data_set = $ARGV[1];}
else {$data_set = "all";}
if($#ARGV >= 2) {$param_set = $ARGV[2];}
else {$param_set = "ideal";}
$param_file = "parameters_theta_${sim_name}_${fit_detector}_${param_set}.dat";
print "Using fit of theta for detector $fit_detector\n";
print "Using parameter file name: $param_file\n";

$input_file = "${sim_name}_rays_${data_set}_${fit_detector}.dat";
print "Input data file name: $input_file\n";

$npar_theta = 7;
@param_label_theta = qw(b_t0 b_a0 b_a1 b_a2 b_b0 b_b1 b_b2);
$npar_phi = 6;
@param_label_phi = qw(c_0 c_1 c_2 c_3 c_4 c_5);

$print_parameters = 1;

# First read in the data
open( IN, $input_file) || die "Cannot open file $input_file\n";
$n = -1;
while(<IN>)
	{
	next if /^#/;
	chomp;
	@data = split(/\s+/);
	if($data[0] =~ /Event:/)
		{ $n++; $ev[$n] = $data[1];} # a new event
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

$d = $fit_detector;
	$param_file_theta = "parameters_theta_geant_${d}_ideal.dat";
	# read in the parameter file for this detector
	open( IN, $param_file_theta) || die "Cannot open file $param_file_theta\n";
	for($i = 0; $i < $npar_theta; $i++)
		{ $par_theta[$d][$i] = 0;}
	
	while(<IN>)
		{
		chomp;
		@data = split(/\s+/);
		for($i = 0; $i < $npar_theta; $i++)
			{
			if($data[0] eq $param_label_theta[$i])
				{
				$par_theta[$d][$i] = $data[1];
				}
			}
		}
	close(IN);
	if($print_parameters)
		{
		print "Theta parameters for detector $d ....\n";
		for($i = 0; $i < $npar_theta; $i++)
			{
			printf "  $param_label_theta[$i]  =  %.4g\n", $par_theta[$d][$i];
			}
		}
	$param_file_phi = "parameters_phi_geant_${d}_ideal.dat";
	# read in the parameter file for this detector
	open( IN, $param_file_phi) || die "Cannot open file $param_file_phi\n";
	for($i = 0; $i < $npar_phi; $i++)
		{ $par_phi[$d][$i] = 0;}
	
	while(<IN>)
		{
		chomp;
		@data = split(/\s+/);
		for($i = 0; $i < $npar_phi; $i++)
			{
			if($data[0] eq $param_label_phi[$i])
				{
				$par_phi[$d][$i] = $data[1];
				}
			}
		}
	close(IN);
	if($print_parameters)
		{
		print "Phi parameters for detector $d ....\n";
		for($i = 0; $i < $npar_phi; $i++)
			{
			printf "  $param_label_phi[$i]  =  %.4g\n", $par_phi[$d][$i];
			}
		}

# Write a file with the fitted function
# and make a histogram
$diff_min = -15.; $diff_max = 15., $diff_bin = 0.4;
$h_num = int(($diff_max - $diff_min)/$diff_bin);
for($h =0; $h < $h_num; $h++)
	{
	$hist[$h] = 0;
	$hist_diff[$h] = $diff_min + $diff_bin*$h;
	}
$outfile = "reconstruct_polar_${sim_name}_${data_set}_${fit_detector}_${param_set}.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";

$sum_x = 0; $sum_2x = 0;
$sum_x_h = 0; $sum_2x_h = 0;
$num_sum_h = 0;
for($n = 0; $n < $num; $n++)
	{
	$theta = function_theta($d, $x_out[$n], $th_out[$n]);
	$phi = function_phi($d, $x_out[$n], $y_out[$n], $th_out[$n], $ph_out[$n]);
	$x = tan($theta/1000.);
	$y = tan($phi/1000.);
	$r = sqrt($x**2 + $y**2);
	$polar_fit = atan($r) * 180. / 3.14159265;
	$x = tan($th_in[$n]/1000.);
	$y = tan($ph_in[$n]/1000.);
	$r = sqrt($x**2 + $y**2);
	$polar_in = atan($r) * 180. / 3.14159265;
	
	print OUT "$polar_in $polar_fit\n";
	$polar_diff = $polar_fit - $polar_in;
	$h = int(($polar_diff - $diff_min)/$diff_bin);
	if($h >= 0 && $h < $h_num)
		{
		$hist[$h]++;
		$sum_x_h += $polar_diff;
		$sum_2x_h += $polar_diff**2;
		$num_sum_h++;
		}
	$sum_x += $polar_diff;
	$sum_2x += $polar_diff**2;
	}
close(OUT);
$outfile = "histogram_polar_${sim_name}_${data_set}_${fit_detector}_${param_set}.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";
for($h =0; $h < $h_num; $h++)
	{
	print OUT "$hist_diff[$h] $hist[$h]\n";
	}
close(OUT);
$average_diff = $sum_x/$num;
$std_diff = $sum_2x - ($sum_x**2)/$num;
$std_diff = sqrt($std_diff/($num-1));
print "Fit to $sim_name data set $data_set\nDetector $fit_detector:\n";
printf "Average difference = %.2f deg Std Dev = %.2f deg\n",$average_diff,$std_diff;
$average_diff = $sum_x_h/$num_sum_h;
$std_diff = $sum_2x_h - ($sum_x_h**2)/$num_sum_h;
$std_diff = sqrt($std_diff/($num_sum_h-1));
printf "In range $diff_min to $diff_max degrees:\n";
printf "Average difference = %.2f deg Std Dev = %.2f deg\n",$average_diff,$std_diff;


exit;

sub function_theta
	{
	my ($val, $x, $p);
	my ($B1, $B2, $B3, $B4, $B5, $B6, $B7);
	$d = $_[0]; $x = $_[1]; $p = $_[2]; 

	$B1 = $par_theta[$d][0] - $par_theta[$d][1]*$par_theta[$d][4];
	$B2 = $par_theta[$d][1];
	$B3 = $par_theta[$d][2];
	$B4 = $par_theta[$d][3];
	$B5 = -($par_theta[$d][2]*$par_theta[$d][4] + $par_theta[$d][1]*$par_theta[$d][5]);
	$B6 = -($par_theta[$d][3]*$par_theta[$d][4] + $par_theta[$d][2]*$par_theta[$d][5] + $par_theta[$d][1]*$par_theta[$d][6]);
	$B7 = -($par_theta[$d][3]*$par_theta[$d][5] + $par_theta[$d][2]*$par_theta[$d][6]);
	$B8 = -$par_theta[$d][3]*$par_theta[$d][6];

	$val = $B1 +$B2*$p +$B3*$p*$x +$B4*$p*$x**2
		+$B5*$x +$B6*$x**2 +$B7*$x**3 +$B8*$x**4;
	return $val;
	}
sub function_phi
	{
	my ($val, $x, $y, $t, $p);
	$d = $_[0];
	$x = $_[1]; $y = $_[2];
	$t = $_[3]; $p = $_[4];

	$val = $par_phi[$d][0]
		+ $par_phi[$d][1]*$p
		+ $par_phi[$d][2]*$p**2
		+ $par_phi[$d][3]*$t
		+ $par_phi[$d][4]*$x
		+ $par_phi[$d][5]*$y;
	return $val;
	}

