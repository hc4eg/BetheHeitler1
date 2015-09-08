#!/usr/bin/perl
use Math::Trig;

$sim_name = "geant";
if($#ARGV >= 0) {$fit_detector = $ARGV[0];}
else {$fit_detector = 1;}
if($#ARGV >= 1) {$data_set = $ARGV[1];}
else {$data_set = "all";}

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
	elsif ($data[0] =~ /Monitor:/)
		{
		if($data[1] =~ /e\+/) {$ip = 0;}
		elsif($data[1] =~ /e-/) {$ip = 1;}
		else {next;}
		$hit_m[$ip][$n] = 1;
		$energy_m[$ip][$n] = $data[2];
		$x_m[$ip][$n] = $data[3];
		$y_m[$ip][$n] = $data[4];
		$th_m[$ip][$n] = $data[5];
		$ph_m[$ip][$n] = $data[6];
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

# Write a file with the fitted function
# and make a histogram
$diff_min = -15.; $diff_max = 15., $diff_bin = 0.4;
$h_num = int(($diff_max - $diff_min)/$diff_bin);
for($h =0; $h < $h_num; $h++)
	{
	$hist[$h] = 0;
	$hist_diff[$h] = $diff_min + $diff_bin*$h;
	}
$outfile = "reconstruct_monitor_polar_${sim_name}_${data_set}_${fit_detector}.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";

$sum_x = 0; $sum_2x = 0;
$sum_x_h = 0; $sum_2x_h = 0;
$num_sum_h = 0;
for($n = 0; $n < $num; $n++)
	{
	$x = tan($th_m[$d][$n]/1000.);
	$y = tan($ph_m[$d][$n]/1000.);
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
$outfile = "histogram_monitor_polar_${sim_name}_${data_set}_${fit_detector}.dat";
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

