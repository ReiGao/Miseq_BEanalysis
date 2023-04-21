#########################################################################
#	File Name: BaseCall.pl
#	> Author: QiangGao
#	> Mail: qgao@qi-biodesign.com
#	Created Time: Sun 15 Oct 2017 04:41:41 PM CST
#	Modify Time: Tue Sep 20 15:08:38 CST 2022
#########################################################################

#!/usr/bin/perl -w
#use strict;
my $QUALCUTOFF=10;
my ($R1,$out)=@ARGV;
my $name;
if($R1=~/extendedFrags/){
 ($name)=$R1=~/.*\/(.*?)\.extendedFrags/;
}else{
  ($name)=$R1=~/.*\/(.*?)\/R/;
}
if($out){
	$name=$name;
}
my $database="./target.txt";
my %database;
my $outdir="splitfile";
my $cc=`mkdir -p $outdir`;
$gene=uc($gene);
open(IN,"$database") or die "$database is not exists\n";
while(<IN>){
	chomp $_;
	next if($_=~/^#/);
	my ($bar,$tar)=split(/\s+/,$_);
	#print "$bar\t$tar\n";
	$database{$bar}{TARGET}=$tar;
	my $R1=$bar."R1";
	my $fa=$bar."FA";
    $hashout{$R1}=$R1;
    $hashout{$fa}=$fa;
}
close IN;
my $R2=$R1;
$R2=~s/R1\.fastq/R2\.fastq/ if($R2=~/R1.fastq/);
if($R1=~/gz$/){
	open(IN,"gzip -dc $R1|");
}else{
	open(IN,"$R1");
}
my $count=0;
my $passqual=10;
my %seq;
my %hashSEQ;
while($a=<IN>){
	$a=~s/\s+//g;
	my $seq=<IN>;
	my $b=<IN>;
	my $qual=<IN>;
	my ($flag,$seq)=qual($seq,$qual);
	if ($flag==1){
		$passqual+=1;
		$seq{$a}=$seq;
		$hashSEQ{$a}=$a."\n".$seq."\n".$b.$qual;
	}
	$count+=1;
}
close IN;

my @name=sort keys %seq;
my $big=@name;
my $in=0;
my $del=0;
my $edit=0;
my $editC=0;
my $onlyC=0;
my $onlyG=0;
my $editG=0;
my $editCG=0;
my $middC=0;
my $noedit=0;
my $all=0;
my %change;
my %site;
my %site2;
my %all;
my @base=("A","T","C","G");
my %indel;
my %indelT;
my $geneCount;
my ($type1,$type2,$type3,$typ4);
open(OUT,">$out.TTT.fa");
#open(OUT2,">$out.nodatabase.fa");
foreach(@name){
	my $seq=$seq{$_};
	my ($gene)=$seq=~/AGCTCT(\w{6})ACTAGT/;
	if(!$gene){
		next;
	}
	if(!(exists $database{$gene}{TARGET})){
		next;
	}
	$geneCount+=1;
	my ($ss)=$seq=~/TGGGAATG(\w{3})TGCCTT/;
	next if(!$ss);
	my $lenss=length($ss);
	my $tar='NNN';
	my ($tar)=$seq=~/TGGGAATG(\w{3})TGCCTT/;
	$type1+=1;
	next if($tar=~/N/);
	my @tar;
	$tar[0]=$tar;
	$all+=1;
	my $OUTR1=$gene."R1";
	my $FA=$gene."FA";
	if(length $tar[0]==length $database{$gene}{TARGET}){
		$all{$tar[0]}+=1;
		if($tar[0]=~/$database{$gene}{TARGET}/){
			$noedit+=1;
		}else{
			my ($e,$ec,$eg,$miC)=compare($tar[0],$database{$gene}{TARGET});
			$middC+=$miC;
			print OUT "Seq\t$ss\n"if($e>0);
			print OUT "Get\t$gene\t$database{$gene}{TARGET}\t$tar[0]\t$e\n" if($e>0);
			$edit+=1 if($e>0);
			$editC+=1 if($ec>0);
			$onlyC+=1 if($ec>0 and $eg==0);
			$editG+=1 if($eg>0);
			$onlyG+=1 if($eg>0 and $ec==0);
			$editCG+=1 if($eg>0 and $ec>0);
		}
	}elsif(length $tar[0] > length $database{$gene}{TARGET}){
		$indel{$tar[0]}+=1;
		$indelT{$tar[0]}="INS";
		print "INS=$tar[0]\n";
		$in+=1;
	}elsif(length $tar[0]<length $database{$gene}{TARGET}){
		$indel{$tar[0]}+=1;
		$indelT{$tar[0]}="DEL";
		print "seq=$seq\n$gene=$database{$gene}{TARGET}) DEL=$tar[0]\n";
		$del+=1;
	}else{
		print "ERRO=$tar[0]\n";
	}
}
if($all==0){
	print "ERR:no LR find $all=0\n";
	print "$out\t$count\t$passqual\t$all\n";
	exit;
}
my $efr=int($all/$count*10000)/100;
my $gfr=int($geneCount/$count*10000)/100;
my $edr=int($edit/$all*10000)/100;
my $edcr=int($editC/$all*10000)/100;
my $oncr=int($onlyC/$all*10000)/100;
my $ongr=int($onlyG/$all*10000)/100;
my $edcgr=int($editCG/$all*10000)/100;
my $edgr=int($editG/$all*10000)/100;
my $midCr=int($middC/$all*10000)/100;
my $type1r=int($type1/$count*10000)/100;
my $idr=int(($in+$del)/$all*10000)/100;
print "Target	AllReads	PassQC	FindLR	FindLRuse	TotalMatchP	TotalMatchP(%)	找到定位边界的Reads	有效比例(%)	没有发生任何编辑的序列数目	单碱基编辑的序列	C>T	C>T(%)	G>A	G>A(%)	bothCG	bothCG(%)	onlyC	onlyC(%)	onlyG	onlyG(%)	插入编辑序列数目	缺失编辑的序列数目	单碱基编辑率(%)	插入缺失编辑率(%)\tMidC\tMidC(%)\n";
print "$out\t$count\t$passqual\t$geneCount\t$gfr\t$type1\t$type1r\t$all\t$efr\t$noedit\t$edit\t$editC\t$edcr\t$editG\t$edgr\t$editCG\t$edcgr\t$onlyC\t$oncr\t$onlyG\t$ongr\t$in\t$del\t$edr\t$idr\t$middC\t$midCr\n";
close OUT;

