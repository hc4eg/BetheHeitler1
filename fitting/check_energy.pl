#!/usr/bin/perl
$sim_name = "geant";
$npar = 10;
@param_label = qw(a_0 a_x a_2x a_3x a_th a_2th a_3th a_xth a_2xth a_x2th);
$central_energy = 30.;

if($#ARGV >= 0) {$fit_detector = $ARGV[0];}
else {$fit_detector = 1;}
if($#ARGV >= 1) {$data_set = $ARGV[1];}
else {$data_set = "all";}
if($#ARGV >= 2) {$param_set = $ARGV[2];}
else {$param_set = "";}
$param_file = "parameters_energy_${sim_name}_${fit_detector}_${param_set}.dat";
print "Using fit of energy for detector $fit_detector\n";
print "Using parameter file name: $param_file\n";

$input_file = "${sim_name}_rays_${data_set}_${fit_detector}.dat";
print "Input data file name: $input_file\n";

## First read in the data
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
$diff_min = -5.; $diff_max = 5., $diff_bin = 0.2;
$h_num = int(($diff_max - $diff_min)/$diff_bin);
for($h =0; $h < $h_num; $h++)
	{ $hist[$h] = 0;
	$hist_diff[$h] = $diff_min + $diff_bin*$h;}
$outfile = "reconstruct_energy_${sim_name}_${data_set}_${fit_detector}_${param_set}.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";

$sum_x = 0; $sum_2x = 0;
$sum_x_h = 0; $sum_2x_h = 0;
$num_sum_h = 0;
for($i = 0; $i < $num; $i++)
	{
	$fit = function($x_out[$i], $th_out[$i] );
	$eng_fit = $central_energy * (1. + $fit/100.);
	print OUT "$energy_in[$i] $eng_fit\n";
	$eng_diff = $eng_fit - $energy_in[$i];
	$h = int(($eng_diff - $diff_min)/$diff_bin);
	if($h >= 0 && $h < $h_num)
		{
		$hist[$h]++;
		$sum_x_h += $eng_diff;
		$sum_2x_h += $eng_diff**2;
		$num_sum_h++;
		}
	$sum_x += $eng_diff;
	$sum_2x += $eng_diff**2;
	}
close(OUT);
$outfile = "histogram_energy_${sim_name}_${data_set}_${fit_detector}_${param_set}.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";
for($h =0; $h < $h_num; $h++)
	{
	print OUT "$hist_diff[$h] $hist[$h]\n";
	}
close(OUT);
$average_diff = $sum_x/$num;
$std_diff = $sum_2x - ($sum_x**2)/$num;
$std_diff = sqrt($std_diff/($num-1));
printf "Average difference = %.2f MeV Std Dev = %.2f MeV\n",$average_diff,$std_diff;
$average_diff = $sum_x_h/$num_sum_h;
$std_diff = $sum_2x_h - ($sum_x_h**2)/$num_sum_h;
$std_diff = sqrt($std_diff/($num_sum_h-1));
printf "In range $diff_min to $diff_max MeV:\n";
printf "Average difference = %.2f MeV Std Dev = %.2f MeV\n",$average_diff,$std_diff;
exit;

sub function
	{
	my ($val, $x, $p);
	$x = $_[0]; $p = $_[1]; 
	$val = $par[0] + $par[1]*$x  + $par[2]*$x*$x + $par[3]*$x**3
		+ $par[4]*$p + $par[5]*$p*$p + $par[6]*$p**3
		+ $par[7]*$x*$p + $par[8]*$x*$x*$p + $par[9]*$x*$p*$p;
	return $val;
	}
