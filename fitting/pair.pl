#!/usr/bin/perl
#
use Math::Trig;
$central_energy = 30.;
$input_file = "geant_rays_photons.dat";
$print_parameters = 1;

$npar_energy = 10; 
@param_label_energy = qw(a_0 a_x a_2x a_3x a_th a_2th a_3th a_xth a_2xth a_x2th);
$npar_theta = 7;
@param_label_theta = qw(b_t0 b_a0 b_a1 b_a2 b_b0 b_b1 b_b2);
$npar_phi = 6;
@param_label_phi = qw(c_0 c_1 c_2 c_3 c_4 c_5);


for($d = 0; $d < 2; $d++)
	{
	$param_file_energy = "parameters_energy_geant_${d}_ideal.dat";
	# read in the parameter file for this detector
	open( IN, $param_file_energy) || die "Cannot open file $param_file_energy\n";
	for($i = 0; $i < $npar_energy; $i++)
		{ $par_energy[$d][$i] = 0;}
	
	while(<IN>)
		{
		chomp;
		@data = split(/\s+/);
		for($i = 0; $i < $npar_energy; $i++)
			{
			if($data[0] eq $param_label_energy[$i])
				{
				$par_energy[$d][$i] = $data[1];
				}
			}
		}
	close(IN);
	if($print_parameters)
		{
		print "Energy parameters for detector $d ....\n";
		for($i = 0; $i < $npar_energy; $i++)
			{
			printf "  $param_label_energy[$i]  =  %.4g\n", $par_energy[$d][$i];
			}
		}
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
	}

# now read the event-by-event data
open( IN, $input_file) || die "Cannot open file $input_file\n";
$n = -1;
while(<IN>)
	{
	next if /^#/;
	chomp;
	@data = split(/\s+/);
	if($data[0] =~ /Event:/) # a new event
		{
		$n++;
		$ev[$n] = $data[1];
		$hit[0][$n] = $hit[1][$n] = 0;
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
		$d = $data[1];
		$hit[$d][$n] = 1;
		}
	elsif ($data[0] =~ /VDC:/)
		{
		$x_out[$d][$n] = $data[1];
		$y_out[$d][$n] = $data[2];
		$th_out[$d][$n] = $data[3];
		$ph_out[$d][$n] = $data[4];
		}
	}
close(IN);

$num = $n + 1;
print "Found $num pair production events\n";

$polar_min = -10.; $polar_max = 10.; $polar_bin = 0.2;
$h_num_polar = int(($polar_max - $polar_min)/$polar_bin);
for($h =0; $h < $h_num_polar; $h++)
	{ $hist_polar[$h] = 0;
	$hist_polar_angle[$h] = $polar_min + $polar_bin*$h;}

$diff_min = -50.; $diff_max = 50., $diff_bin = 2.;
$h_num_diff = int(($diff_max - $diff_min)/$diff_bin);
for($h =0; $h < $h_num_diff; $h++)
	{ $hist_diff[$h] = 0;
	$hist_diff_eng[$h] = $diff_min + $diff_bin*$h;}

$sum_min = 40.; $sum_max = 80., $sum_bin = 0.5;
$h_num_sum = int(($sum_max - $sum_min)/$sum_bin);
for($h =0; $h < $h_num_sum; $h++)
	{ $hist_sum[$h] = 0;
	$hist_sum_eng[$h] = $sum_min + $sum_bin*$h;}

$open_min = 0.; $open_max = 20., $open_bin = 0.5;
$h_num_open = int(($open_max - $open_min)/$open_bin);
for($h =0; $h < $h_num_open; $h++)
	{ $hist_open[$h] = 0;
	$hist_open_angle[$h] = $open_min + $open_bin*$h;}

$a_open_min = 2.5; $a_open_max = 27.5; $a_open_bin = 5.;
$a_diff_min = -30.; $a_diff_max = 30.; $a_diff_bin = 4.;
$a_num_open = int(($a_open_max - $a_open_min)/$a_open_bin);
$a_num_diff = int(($a_diff_max - $a_diff_min)/$a_diff_bin);
for($h =0; $h < $a_num_open; $h++)
	{$a_open_angle[$h] = $a_open_min + $a_open_bin*($h + 0.5);}
for($h =0; $h < $a_num_diff; $h++)
	{$a_eng_diff[$h] = $a_diff_min + $a_diff_bin*($h + 0.5);}
for($i = 0; $i < $a_num_open; $i++)
	{ for($j = 0; $j < $a_num_diff; $j++)
		{ $a_n[$i][$j] = 0;} }

$sig_0 = 2.51; $sig_0 = 2.54; #from Uranium target
$angle_tol = 2.* sqrt( $sig_0**2 + $sig_1**2); #degrees
printf "Tolerance angle = %.2f\n", $angle_tol;

# loop over all events looking for a pair production event
$npair = 0;
for($n = 0; $n < $num; $n++)
	{
	if($hit[0][$n] && $hit[1][$n])
		{
		$npar++;
		for($d = 0; $d < 2; $d++)
			{
			$fit = function_energy($d, $x_out[$d][$n], $th_out[$d][$n]);
			$eng[$d] = $central_energy * (1. + $fit/100.);
			$theta[$d] = function_theta($d, $x_out[$d][$n], $th_out[$d][$n]);
			$phi[$d] = function_phi($d, $x_out[$d][$n], $y_out[$d][$n], $th_out[$d][$n], $ph_out[$d][$n]);
			$x[$d] = tan($theta[$d]/1000.);
			$y[$d] = tan($phi[$d]/1000.);
			$r[$d] = sqrt($x[$d]**2 + $y[$d]**2 + 1.);
			$an[$d] = acos(1./$r[$d]) * 180. / 3.14159265;
			}
		$energy_sum = $eng[0] + $eng[1];
		$energy_diff = $eng[1] - $eng[0]; # electron - positron
		#print "Event $ev[$n]: eng[0] = $eng[0] eng[1] = $eng[1]\n";
		$h = int(($energy_diff - $diff_min)/$diff_bin);
		if($h >= 0 && $h < $h_num_diff)
			{
			$hist_diff[$h]++;
			$sum_x_diff += $energy_diff;
			$sum_2x_diff += $energy_diff**2;
			$num_diff++;
			}
		$h = int(($energy_sum - $sum_min)/$sum_bin);
		if($h >= 0 && $h < $h_num_sum)
			{
			$hist_sum[$h]++;
			$sum_x_sum += $energy_sum;
			$sum_2x_sum += $energy_sum**2;
			$num_sum++;
			}
		$dot = $x[0]*$x[1] + $y[0]*$y[1] + 1.;
		$ang = acos( $dot/($r[0]*$r[1]) );
		$angle = $ang * 180. / 3.14159265;
		
		$h = int(($angle - $open_min)/$open_bin);
		if($h >= 0 && $h < $h_num_open)
			{
			$hist_open[$h]++;
			$sum_x_open += $angle;
			$sum_2x_open += $angle**2;
			$num_open++;
			}

		# asymmetry stuff
		if(abs($an[0] - $an[1]) < $angle_tol)
			{
			$i = int(($an[1] - $a_open_min)/$a_open_bin);
			if($i >= 0 && $i < $a_num_open)
				{
				$j = int(($energy_diff - $a_diff_min)/$a_diff_bin);
				if($j >= 0 && $i < $a_num_diff)
					{
					$a_n[$i][$j]++;
					$npair++;
					}
				}
			}
		}
	}
$average_diff = $sum_x_diff/$num_diff;
$std_diff = $sum_2x_diff - ($sum_x_diff**2)/$num_diff;
$std_diff = sqrt($std_diff/($num_diff-1));
printf "Average difference = %.2f MeV Std Dev = %.2f MeV\n",$average_diff,$std_diff;
$average_sum = $sum_x_sum/$num_sum;
$std_sum = $sum_2x_sum - ($sum_x_sum**2)/$num_sum;
$std_sum = sqrt($std_sum/($num_sum-1));
printf "Average sum = %.2f MeV Std Dev = %.2f MeV\n",$average_sum,$std_sum;
$average_open = $sum_x_open/$num_open;
$std_open = $sum_2x_open - ($sum_x_open**2)/$num_open;
$std_open = sqrt($std_open/($num_open-1));
printf "Average Opening angle = %.2f deg Std Dev = %.2f deg\n",$average_open,$std_open;

print "Events surviving symmetry requirement = $npair\n";

$outfile = "pair_energy_sum.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";
for($h =0; $h < $h_num_sum; $h++)
	{
	print OUT "$hist_sum_eng[$h] $hist_sum[$h]\n";
	}
close(OUT);
$outfile = "pair_energy_difference.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";
for($h =0; $h < $h_num_diff; $h++)
	{
	print OUT "$hist_diff_eng[$h] $hist_diff[$h]\n";
	}
close(OUT);
$outfile = "pair_open_angle.dat";
open(OUT, ">$outfile") || die "Cannot open file $outfile\n";
for($h =0; $h < $h_num_open; $h++)
	{
	print OUT "$hist_open_angle[$h] $hist_open[$h]\n";
	}
close(OUT);

# calculate asymmetries
for($i = 0; $i < $a_num_open; $i++)
	{
	for($j = 0; $j < $a_num_diff; $j++)
		{
		$jj = $a_num_diff - $j - 1;
		$sp = $a_n[$i][$j];
		$sm = $a_n[$i][$jj];
		$s = ($sp + $sm); 
		if($sp > 0. && $sm > 0.)
			{
			$asym[$i][$j] = ($sp - $sm)/$s; 
			$asym_err[$i][$j] = 2.*sqrt($sp*$sm*$s)/$s/$s;
			}
		else
			{
			$asym[$i][$j] = 0.;
			$asym_err[$i][$j] = 0.;
			}
		}
	$outfile = "pair_asymmetry_$a_open_angle[$i].dat";
	open(OUT, ">$outfile") || die "Cannot open file $outfile\n";	
	for($j = 0; $j < $a_num_diff; $j++)
		{
		if($asym_err[$i][$j] > 0.) {
		printf OUT "%.2f %.4f %.4f\n", $a_eng_diff[$j], $asym[$i][$j], $asym_err[$i][$j]; }
		}
	close(OUT);
	}



sub function_energy
	{
	my ($val, $x, $p, $d);
	$d = $_[0]; $x = $_[1]; $p = $_[2]; 
	$val = $par_energy[$d][0]
		+ $par_energy[$d][1]*$x
		+ $par_energy[$d][2]*$x*$x
		+ $par_energy[$d][3]*$x**3
		+ $par_energy[$d][4]*$p
		+ $par_energy[$d][5]*$p*$p
		+ $par_energy[$d][6]*$p**3
		+ $par_energy[$d][7]*$x*$p
		+ $par_energy[$d][8]*$x*$x*$p
		+ $par_energy[$d][9]*$x*$p*$p;
	return $val;
	}

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

