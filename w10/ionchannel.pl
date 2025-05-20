#!/usr/bin/perl -w
# ionchannel.pl

use strict;
use POSIX;
use Class::Struct;
use lib '/Users/erivas/projects/R-scape/scripts';
#use lib '/Users/rivase/projects/R-scape/scripts';
use FUNCS;
use constant GNUPLOT => '/opt/homebrew/bin/gnuplot';
#use constant GNUPLOT => '/usr/local/bin/gnuplot';
#use constant GNUPLOT => '/opt/local/bin/gnuplot';


use vars qw ($opt_t $opt_N $opt_P $opt_R $opt_T);  # required if strict used
use Getopt::Std;
getopts ('tN:P:R:T:');

# Print a helpful message if the user provides no input file.
if (!@ARGV) {
        print "usage:  ionchannel.pl [options] <ko> <kc> <file> \n\n";
        print "options:\n";
        exit;
}

my $ko     = shift;   # rate of channel open
my $kc     = shift;   # rate of channel closed
my $file   = shift;

my $dogillespie = 0;
my $N = 0;
if ($opt_N) { $N = $opt_N; }

my $T  = 1000;
if ($opt_T) { $T = $opt_T; }
my $KT = 1000;
my $K  = 1000;

my $cfile = "$file.c";
my $sran = rand();
my $PO_nought  = 1.0;   ## P(Open) at t=0
if ($sran < 0.5) { $PO_nought  = 0.0; }

open(CFILE, ">$cfile");
for (my $k = 0; $k < $KT; $k ++) {
    my $time = $k/$KT * $T;
    
    my $PO = channel_continuous($ko, $kc, $time, $PO_nought);
    printf CFILE "$time $PO\n";
    
}
close(CFILE);

my $M = $KT;
my $aveplotfile = "$file.ave";
create_averages($M, $aveplotfile);

my @plotfile;
my @num;
my @Copen;
for (my $n = 0; $n < $N; $n ++) {
    $plotfile[$n] = "$file.$n";
    if ($dogillespie) { gillespie ($T, $KT, $ko, $kc, $PO_nought, \@Copen, \@num, $plotfile[$n]); }
    else              { MonteCarlo($T, $KT, $ko, $kc, $PO_nought, \@Copen, \@num, $plotfile[$n]); }
}
plot_outfile($file, $N, \@plotfile, $aveplotfile);


sub gillespie {
    my ($T, $KT, $ko, $kc, $PO_nought, $Copen_ave_ref, $num_ref, $outfile) = @_;
    
    if ($outfile) { open (OUT, ">$outfile") || die; }
    my $Copen = $PO_nought;

    my $time = 0;
    my $cur_k = ceil($time*$KT/$T);
    my $prv_time = $time;
    my $prv_k    = $cur_k;

    while ($time <= $T) {
	my $c0 = ($Copen == 1)? $kc : $ko;
	
	my $lambda = $c0;
        my $tw     = sample_exponential($lambda);
	$prv_time  = $time;
	$time     += $tw;
	$prv_k     = ceil($prv_time*$KT/$T);
	$cur_k     = ceil(    $time*$KT/$T);

	for (my $k = $prv_k; $k < $cur_k; $k ++) {
	    my $tt = $k/$KT*$T;
	    $num_ref->[$k] ++;
	    $Copen_ave_ref->[$k] += $Copen;
	    if ($outfile) {
		printf OUT "%.2f %f\n", $tt, $Copen;
	    }
	}
	
	my $newC;
	if    ($Copen == 1) { $newC = 0; }
	elsif ($Copen == 0) { $newC = 1; }
	$Copen = $newC;
    }
    
    if ($outfile) { close(OUT); }
}

sub MonteCarlo {
    my ($T, $KT, $ko, $kc, $PO_nought, $Copen_ave_ref, $num_ref, $outfile) = @_;
    
    if ($outfile) { open (OUT, ">$outfile") || die; }
    my $Copen = $PO_nought;
    my $time = 0;
    for (my $k = 0; $k < $KT; $k ++) {
	my $deltat = $T/$KT;
	my $c0 = ($Copen == 1)? $kc : $ko;
	my $sran = rand()/$deltat;

	my $newC;
	if ($sran < $c0) {
	    if    ($Copen == 1) { $newC = 0; }
	    elsif ($Copen == 0) { $newC = 1; }
	    $Copen = $newC;
	}
	
       	#printf     "%.2f %f\n", $time,     $Copen;
	if ($outfile) {
	    printf OUT "%.2f %f\n", $time,         $Copen;
	    printf OUT "%.2f %f\n", $time+$deltat, $Copen;
	}

	$num_ref->[$k] ++;
	$Copen_ave_ref->[$k] += $Copen;
	$time                += $deltat;
    }

    if ($outfile) { close(OUT); }
}

sub create_averages {
    ($M, $aveplotfile) = @_;
    my @Copen_ave;
    my @num;
    for (my $k = 0; $k < $KT; $k ++) {
	$Copen_ave[$k] = 0;
	$num[$k]       = 0;
    }
    for (my $m = 0; $m < $M; $m ++) {
	my $Copen;
	if ($dogillespie) { gillespie ($T, $KT, $ko, $kc, $PO_nought, \@Copen_ave, \@num, ); }
	else              { MonteCarlo($T, $KT, $ko, $kc, $PO_nought, \@Copen_ave, \@num, ); }
    }
    open (AVE, ">$aveplotfile") || die;
    for (my $k = 0; $k < $KT; $k++) {
	if ($num[$k] > 0) { 
	    my $time = $k/$KT*$T;
	    $Copen_ave[$k] /= $num[$k]; 
	    printf AVE "$time $Copen_ave[$k] $num[$k]\n";
	}
    }
    close(AVE);
}


sub channel_continuous {
    my ($kb, $ku, $time, $P1no) = @_;

    my $srate = $kb+$ku;
    my $rrate = ($srate > 0)? $kb/$srate : 0;
    
    my $P1 = $rrate + ($P1no-$rrate) * exp(-$srate*$time);
    return $P1;
}


sub plot_outfile {
    my ($file, $N, $outfile_ref, $aveplotfile) = @_;

    my $seeplots = 1;
   
    my $title = "ko = $ko kc = $kc";
    my $psfile = "$file.ps";
    my $pdffile = $psfile;
    if ($pdffile =~ /^(\S+).ps$/) { $pdffile = "$1.pdf"; }

    open(GP,'|'.GNUPLOT) || die "Gnuplot: $!";
    
    print GP "set terminal postscript color solid 14\n";
    FUNCS::gnuplot_define_styles (*GP);

    print GP "set output '$psfile'\n";
    #print GP "set nokey\n";
    print GP "set xlabel 'time'\n";
    #print GP "set xrange [$xleft:$xright]\n";
    #print GP "set yrange [-0.01:1.01]\n";
    print GP "set style fill solid noborder \n";

    #print GP "set title \"$title\\n\\n$key\"\n";
    print GP "set title '$title'\n";
    
    my $key = "";
    
    my $m = 4;
    $title = "Ion Channel Gating";
    print GP "set ylabel 'probability '\n";
    print GP "set yrange [-0.1:1.1]\n";
    
    my $cmd  =  "'$cfile'        using 1:2 with lines title '' ls $m, "; $m ++;
    for (my $n = 0; $n < $N-1; $n ++) {
	if ($m == 5) { $m ++; }
	$cmd  .= "'$outfile_ref->[$n]'  using 1:2 with lines title '' ls $m, "; $m ++;
    }
    if ($N > 0) { 
	if ($m == 5) { $m ++; }
	if ($aveplotfile) { $cmd  .= "'$outfile_ref->[$N-1]'  using 1:2 with lines title '' ls $m, "; }
	else              { $cmd  .= "'$outfile_ref->[$N-1]'  using 1:2 with lines title '' ls $m";   }
    }
    if ($aveplotfile) {
	if ($N > 0) { $cmd    .=  "'$aveplotfile'  using 1:2 with lines title 'averages' ls 5, "; }
	else        { $cmd    .=  "'$aveplotfile'  using 1:2 with lines title 'averages' ls 5, "; }
    }
    $cmd  .=  "'$cfile'        using 1:2 with lines title 'deterministic' ls 4 ";
    print    "plot $cmd\n";
    print GP "plot $cmd\n";

    close (GP);
    
    system ("ps2pdf $psfile\n"); 
    system("rm $psfile\n");
    if ($seeplots) { system ("open $pdffile&\n"); }

}

# lambda exp(-lambda*t)
#
sub sample_exponential {
    my ($lambda) = @_;

    if ($lambda == 0) { return 0; }
    
    my $sran = rand();
    my $tw = -1.0/$lambda * log($sran);

    return $tw;
}
