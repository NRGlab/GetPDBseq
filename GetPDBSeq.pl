#!/usr/bin/perl -w

$usage = $0." -[pl] <pdb or list file name> -c <chain> -a <atom seq> -b <outname> -x <create index>";
$list="";
$atom_flag=0;
$base="";
$chn="";
$create_index=0;
$error = 0;
$step=3;

@pdbs=();
while(@ARGV){
    $arg=shift @ARGV;
    #print $arg;<STDIN>;
    if($arg eq "-p"){push(@pdbs,lc(shift @ARGV));}
    if($arg eq "-l"){$list=shift @ARGV;}
    if($arg eq "-c"){$chn=shift @ARGV;}
    if($arg eq "-a"){$atom_flag=1;}
    if($arg eq "-b"){$base=shift @ARGV;}
    if($arg eq "-x"){$create_index=1;}
    if($arg eq "-h"){die "$usage\n";}
    if($arg eq "-n"){$step = 1;}
}
if($base ne ""){$base .= "-";} 
if($list ne ""){
    open(IN,$list) || die "cannot open PDB list $list for read in GetPDBSeqs\n";
    while($line = <IN>){
	chomp $line;
	push(@pdbs,lc($line));
    }
}
close(IN);

if($#pdbs == -1){$error = 1;}
if($error == 1){die "$usage\n";}

foreach $pdb (@pdbs){
    ($id = $pdb) =~ s/\.pdb//;
    #print $id;<STDIN>;
    @tmp=split(/\//,$id);
    #print $tmp[$#tmp];<STDIN>;
    $base{$pdb}=$tmp[$#tmp];
}

if($atom_flag == 0){
    foreach $pdb (@pdbs){print &Seq_from_SEQRES($base{$pdb},$pdb,$chn);}
}else{
    if($create_index==0){
	foreach $pdb (@pdbs){print &Seq_from_ATOM($base{$pdb},$pdb,$chn);}
    }else{
	foreach $pdb (@pdbs){print &Index_from_ATOM($base{$pdb},$pdb,$chn);}
    }
}

#1234567890123456789012345678901234567890123456789012345678901234567890123456#
#         1         2         3         4         5         6         7      #
#1234567890123456789012345678901234567890123456789012345678901234567890123456#
sub Index_from_ATOM(){
    use strict;
    my $base=shift;
    my $pdb = shift;
    my $chn_in = shift;
    my $line;
    my ($res,$chn,$num,$k);
    my $id;
    my %sek;
    my %ids;
    local *IN;
    my $allseqs="";
    my $thisseq;

    if($base eq ""){$base = substr($pdb,0,4)."-";}
    if($pdb !~ /\.pdb$/){$pdb .= ".pdb";}
    open(IN,$pdb) || die "cannot read file $pdb\n";
    
    while($line=<IN>){
	if($line =~ /^ATOM /){
	    $res=substr($line,17,3);
	    ($chn=substr($line,21,1))=~s/\s/\_/;
	    ($num=substr($line,22,4))=~s/\s//g;
	    $id=$res.$chn.$num;
	    if(!exists($ids{$id})){
		push(@{$sek{$chn}},$res.$num);
		$ids{$id}=1;
	    }
	}
    }
    close(IN);
    
    foreach $chn (sort keys %sek){
	if($chn_in ne "" && $chn ne $chn_in){next;}
	$allseqs .= ">$base$chn\n";
	$thisseq="";
	while(@{$sek{$chn}}){
	    $k=shift @{$sek{$chn}};
	    $res=substr($k,0,3);
	    $num=substr($k,3);
	    #print $res,$chn,$num;<STDIN>;
	    $thisseq .= &Translate2One($res).$num.":";	    
	}
	chop $thisseq;
	$allseqs .= $thisseq."\n";
    }
    return $allseqs;
}
#-----------------------------------------------------------------------------
sub Seq_from_ATOM(){
    use strict;
    my $base=shift;
    my $pdb = shift;
    my $chn_in = shift;
    my $line;
    my ($res,$chn,$num,$k);
    my $id;
    my %sek;
    my %ids;
    local *IN;
    my $allseqs="";
    my $thisseq;

    if($base eq ""){$base = substr($pdb,0,4)."-";}
    if($pdb !~ /\.pdb$/){$pdb .= ".pdb";}
    open(IN,$pdb) || die "cannot read file $pdb\n";
    
    while($line=<IN>){
	if($line =~ /^ATOM /){
	    $res=substr($line,17,3);
	    ($chn=substr($line,21,1))=~s/\s/\_/;
	    ($num=substr($line,22,4))=~s/\s//g;
	    $id=$res.$chn.$num;
	    if(!exists($ids{$id})){
		push(@{$sek{$chn}},$res.$num);
		$ids{$id}=1;
	    }
	}
    }
    close(IN);
    
    foreach $chn (sort keys %sek){
	if($chn_in ne "" && $chn ne $chn_in){next;}
	$allseqs .= "> $base$chn\n";
	$thisseq="";
	while(@{$sek{$chn}}){
	    $k=shift @{$sek{$chn}};
	    $res=substr($k,0,3);
	    $num=substr($k,3);
	    #print $res,$chn,$num;<STDIN>;
	    $thisseq .= $res;	    
	}
	$allseqs .= Nice_seq(&Translate2One($thisseq));
    }
    return $allseqs;
}
#-----------------------------------------------------------------------------
sub Seq_from_SEQRES(){
    use strict;

    my $base=shift;
    my $pdb = shift;
    my $chn_in = shift;
    my ($line,$len,$maxlen,$chn,$k);
    my @tmp;
    my %seq=();
    my $file;
    local *IN;

    if($base eq ""){$base = substr($pdb,0,4)."-";}
    if($pdb !~ /\.pdb$/){$pdb .= ".pdb";}
    open (IN,$pdb) || die "cannot open file $pdb for read in sub Seq_from_SEQRES\n";
    while($line=<IN>){
        chomp $line;
        if($line !~ /^SEQRES/){next;}
        $len=70;
        if(length($line) >= $len){
            $maxlen=$len-19;
        }else{
            $maxlen=length($line)-19;
        }
        $chn = substr($line,11,1);
        $seq{$chn} .= substr($line,19,$maxlen)."\n"; 
    }
    close(IN);
    my $allseqs="";
    foreach $k (sort keys %seq){
      #print "**\n",$seq{$k};<STDIN>;
      $chn=$k;
      if($chn_in ne "" && $chn ne $chn_in){next;}
      if($chn eq " "){$chn="_";}
      $line="> $base$chn\n";
      $allseqs .= $line.Nice_seq(&Translate2One($seq{$k}));
    }

    return $allseqs;
}

#----------------------------------------------------------------------------
sub Translate2One(){
    use strict;
    my $seq_ori=shift;
    my $seq_new;
    my $i;
    my $imax;
    our ($step);

    #print "$seq_ori";<STDIN>;
    $seq_ori =~ tr/A-Z//cd;
    #print"step=$step\n";
    if($step == 3){
      while(length($seq_ori) != int(length($seq_ori)/3)*3){$seq_ori .= " ";}
      $imax=length($seq_ori);
      for($i=0;$i<$imax;$i += 3){$seq_new .= Three2One(substr($seq_ori,$i,3));}
    }else{
      $seq_new=$seq_ori;
    }

    #print "<",$seq_ori,">";<STDIN>;
    #print "<",$seq_ori,">";<STDIN>;
    #print "len= ",length($seq_ori),"\n";<STDIN>;
    
    #print "<",$seq_new,">";<STDIN>;
    
    return $seq_new;
}
#----------------------------------------------------------------------------
sub Three2One(){
    use strict;
    my $three=lc(shift);
    my $one="x";

    #print "3=<$three>";
    if($three eq "gly"){$one="g";}
    if($three eq "ala"){$one="a";}
    if($three eq "val"){$one="v";}
    if($three eq "leu"){$one="l";}
    if($three eq "ile"){$one="i";}
    if($three eq "met"){$one="m";}
    if($three eq "pro"){$one="p";}
    if($three eq "phe"){$one="f";}
    if($three eq "trp"){$one="w";}
    if($three eq "ser"){$one="s";}
    if($three eq "thr"){$one="t";}
    if($three eq "asn"){$one="n";}
    if($three eq "gln"){$one="q";}
    if($three eq "tyr"){$one="y";}
    if($three eq "cys"){$one="c";}
    if($three eq "lys"){$one="k";}
    if($three eq "arg"){$one="r";}
    if($three eq "his"){$one="h";}
    if($three eq "asp"){$one="d";}
    if($three eq "glu"){$one="e";}
    if($three eq "  a"){$one="a";} #    A                    Adenine
    if($three eq "  c"){$one="c";} #    C                    Cytosine
    if($three eq "  g"){$one="g";} #    G                    Guanine
    if($three eq "  t"){$one="t";} #    T                    Thymine
    if($three eq "  u"){$one="u";} #    U                    Uracil
    if($three eq "  r"){$one="r";} #    A or G               puRine
    if($three eq "  y"){$one="y";} #    C or T (U)           pYrimidine
    if($three eq "  m"){$one="m";} #    A or C               aMino
    if($three eq "  k"){$one="k";} #    G or T (U)           Keto
    if($three eq "  s"){$one="s";} #    C or G               Strong (triple ‘3 H’ bonds)
    if($three eq "  w"){$one="w";} #    A or T (U)           Weak (double ‘2 H’ bonds)
    if($three eq "  b"){$one="b";} #    C or G or T (U)      not A
    if($three eq "  d"){$one="d";} #    A or G or T (U)      not C
    if($three eq "  h"){$one="h";} #    A or C or T (U)      not G
    if($three eq "  v"){$one="v";} #    A or C or G          not T (U)
    if($three eq "  n"){$one="n";} #    A or C or G or T (U) aNy nucleotide

    #print "1=<$one>\n";<STDIN>;
    return uc($one);
}
#----------------------------------------------------------------------------
sub Nice_seq(){
  use strict;

  my $seq=shift;
  #print "seq in nice\n",$seq;<STDIN>;
  my $n_col=1;
  my $c_len=50;
  my $header="";
  #my $header="SQPROT   ";

  my $len=length($seq);
  my $blocks=int($len/$c_len);
  my $l;
  my $s .= $header;
  my $colcount=0;
  
  for($l=0;$l<$blocks;$l++){
    $s .= substr($seq,$l*$c_len,$c_len)." ";
    $colcount++;
    if($colcount == $n_col){
      chop $s;
      $s .= "\n".$header;
      $colcount=0;
    }
  }
  if($len-$blocks*$c_len > 0){
      $s .= substr($seq,$l*$c_len)."\n";
  }elsif($len-$blocks*$c_len < 0){
      $s = $seq."\n"; 
  }else{
    chop $s;
    $s .="\n";
  }
  
  return $s;
}
#----------------------------------------------------------------------------
