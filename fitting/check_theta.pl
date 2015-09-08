#!/usr/bin/perl
$sim_name = "geant";
if($#ARGV >= 0) {$fit_detector = $ARGV[0];}
else {$fit_detector = 1;}
if($#ARGV >= 1) {$data_set = $ARGV[1];}
else {$data_set = "all";}
if($#ARGV >= 2) {$param_set = $ARGV[2];}
else {$param_set = "";}
$param_file = "parameters_theta_${sim_name}_${fit_detector}_${param_set}.dat";
print "Using fit of theta for detector $fit_detector\n";
print "Using parameter file name: $param_file\n";

$input_file = "${sim_name}_rays_${data_set}_${fit_detector}.dat";
print "Input data file name: $input_file\n";

$npar = 7;
@param_label = qw(b_t0 b_a0 b_a1 b_a2 b_b0 b_b1 b_b2);

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

# read in the parameter file
$param_file = "parameters_theta_${sim_name}_${fit_detector}_${param_set}.dat";
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
			}
		}
	}
close(IN);

print "Input parameters....\n";
for($i = 0; $i < $npar; $i++)
	{
	printf "  $param_label[$i]  =  %.4g\n", $par[$i];
	}

# Write a file with the fitted function
# and make a histogram
$diff_min = -10.; $diff_max = 10., $diff_bin = 0.2;
$h_num = int(($diff_max - $diff_min)/$diff_bin);
for($h =0; $h < $h_num; $h++)
	{
	$hist[$h] = 0;
	$hist_diff[$h] = $diff_min + $diff_bin*$h;
	}
$outfile = "reconstruct_theta_${sim_name}_${data_set}_${fit_detector}_${param_set}.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";

$sum_x = 0; $sum_2x = 0;
$sum_x_h = 0; $sum_2x_h = 0;
$num_sum_h = 0;
for($i = 0; $i < $num; $i++)
	{
	$theta_fit = function($x_out[$i], $th_out[$i] );
	print OUT "$th_in[$i] $theta_fit\n";
	$theta_diff = $theta_fit - $th_in[$i];
	$diff_deg = $theta_diff*180./3.14159/1000.;
	$h = int(($diff_deg - $diff_min)/$diff_bin);
	if($h >= 0 && $h < $h_num)
		{
		$hist[$h]++;
		$sum_x_h += $theta_diff;
		$sum_2x_h += $theta_diff**2;
		$num_sum_h++;
		}
	$sum_x += $theta_diff;
	$sum_2x += $theta_diff**2;
	}
close(OUT);
$outfile = "histogram_theta_${sim_name}_${data_set}_${fit_detector}_${param_set}.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";
for($h =0; $h < $h_num; $h++)
	{
	print OUT "$hist_diff[$h] $hist[$h]\n";
	}
close(OUT);
$average_diff = $sum_x/$num;
$std_diff = $sum_2x - ($sum_x**2)/$num;
$std_diff = sqrt($std_diff/($num-1));
$average_diff_deg = $average_diff*180./3.14159/1000.;
$std_diff_deg = $std_diff*180./3.14159/1000.;
print "Fit to $sim_name data set $data_set\nDetector $fit_detector:\n";
printf "Average difference = %.2f mrad Std Dev = %.2f mrad\n",$average_diff,$std_diff;
printf "Average difference = %.2f deg Std Dev = %.2f deg\n",$average_diff_deg,$std_diff_deg;
$average_diff = $sum_x_h/$num_sum_h;
$std_diff = $sum_2x_h - ($sum_x_h**2)/$num_sum_h;
$std_diff = sqrt($std_diff/($num_sum_h-1));
$average_diff_deg = $average_diff*180./3.14159/1000.;
$std_diff_deg = $std_diff*180./3.14159/1000.;
printf "In range $diff_min to $diff_max degrees:\n";
printf "Average difference = %.2f mrad Std Dev = %.2f mrad\n",$average_diff,$std_diff;
printf "Average difference = %.2f deg Std Dev = %.2f deg\n",$average_diff_deg,$std_diff_deg;


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

