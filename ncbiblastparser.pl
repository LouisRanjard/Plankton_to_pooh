#!/usr/bin/perl

use strict;
use Bio::SearchIO;
use Bio::SeqIO;
#use Bio::DB::GenBank;
use Bio::DB::EUtilities;
use Bio::DB::Taxonomy;
use List::Util qw[min max];


### Usage information
die "Usage: $0 <BLAST-report-file> <number-of-top-hits> <output-file>\n", if (@ARGV != 3);

my ($infile,$numHits,$outfile) = @ARGV;

print "Parsing the BLAST result ...\n";

my $in = Bio::SearchIO->new(-format => 'blast', -file => $infile);

open (OUT,">$outfile") or die "Cannot open $outfile $!";

open (OUTLIST,">$outfile.list") or die "Cannot open $outfile.list $!";
print OUTLIST "label\tnumOtus\t";

#print OUT "query_name\torganism\n";

my @query_nam = ();
my @access_num = ();
my @number_of_hits = ();
my $i = 0;
my $is = 0;
my $numtotal = 0; # total number of accession numbers


### Store the accession numbers of each blast result
while ( my $result = $in->next_result ) {
   	push(@query_nam, $result->query_name); 		# the name of the query sequence
   	push(@number_of_hits, $result->num_hits);	# the number of hits of the query
    	if ( $result->num_hits == 0 ) { # output "no hits found" if there is no hits
		#print "\tNo hits found\n";
		push(@access_num, "No hits found");
    	#} elsif ( $result->num_hits < $numHits ) {
	#	print "\tInsufficient number of hits in the Blast file provided (found " . $result->num_hits . " while " . $numHits . " are required) \n";
	#	exit();
	} else {
		my $count = 0;
		while (my $hit = $result->next_hit) {
			#print $hit->accession . "\t";
			#my $seq = eval{$db_obj->get_Seq_by_acc('JN5123')};
			push(@access_num, $hit->accession);
			#push @{ $access_num[$i] }, $hit->accession; # 2D array
			$count++;
			# flow control for the number of hits needed
			last if ($count == $numHits);
		}
		$numtotal = $numtotal+$count ;
	}
	$is = $i+1;
	print OUTLIST "Otu$is\t";
	$i++;
	#print OUT "\n";
}

#print $i."\n";

### format a taxonomy.list file for mothur
$i=0;
print OUTLIST "\n0.03\t$is\t";
foreach (@query_nam) {
	print $query_nam[$i] . " " . $number_of_hits[$i] . "\n";
	if ($number_of_hits[$i]==0) {
		print OUTLIST "NoBlast" ; # create a dummy sequence name for sequences that didn't get any Blast results
	} else {
		for my $c (1..$number_of_hits[$i]) {
			$_ =~ s/[;:]/_/g; # replace ";:" with "_" in the query name
			print OUTLIST $_."_".$c;
			if ($c<$number_of_hits[$i]) {print OUTLIST ",";}
		}
	}
	$i++;
	print OUTLIST "\t";
}
close OUTLIST;


#print join(", ", @access_num);
#print scalar @access_num;

print "Total number of accession numbers to download: ".$numtotal."\n";
#print $is;

#my $db = Bio::DB::GenBank->new();
#my $seqio = $db->get_Stream_by_acc(\@access_num);
#my $seqio = $db->get_Seq_by_acc(\@access_num);

### retrieve Genbank format for all accession numbers
#my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
#                                       -db      => 'nucleotide',
#                                       -rettype => 'gb',
#                                       -email   => 'l.ranjard@auckland.ac.nz',
#                                       -term    => join(',',@access_num) );
#                                       -id      => \@access_num); # search accession number in the field ID (ok with EUtilities)
#my $file = 'tmp.gb';
# dump HTTP::Response content to a file (not retained in memory)
#$factory->get_Response(-file => $file);
#my $seqio = Bio::SeqIO->new(-file   => $file,
#                            -format => 'genbank');

### same with breaking down the number of numbers to fetch at a time
#open(my $file, '>>', 'tmps.gb');
#my $file1 = 'tmp1.gb';
#my $retstartnum = 1 ;
#my $retmaxnum = 10 ;
#while ( $retstartnum<$numtotal ) {
#	print $retstartnum . " " . $retmaxnum . "\n";
#	my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
#		                               -db      => 'nucleotide',
#		                               -rettype => 'gb',
#		                               -email   => 'l.ranjard@auckland.ac.nz',
#		                               -id      => \@access_num,
#					       -retstart=> $retstartnum,
#					       -retmax  => $retmaxnum );
#	print $file $factory->get_Response() ;
#	my $file0 = 'tmp'.$retstartnum.'.gb' ;
#	$factory->get_Response(-file => $file0);
#	#cat tmp0.gb >> tmp1.gb ;
#	$retstartnum = $retstartnum+$retmaxnum ;
#}

### same with breaking down the number of sequences to (i) submit and (ii) retrieve, following http://www.bioperl.org/wiki/HOWTO:EUtilities_Cookbook
### max for query (querymax) 10,000
### max for retrieve (retmax) 500
my $querymax = min(1000, $numtotal);
#my $querymax = min(200, $numtotal); # debug
my $retry = 0;
my ($year, $mon, $mday, $hour, $min, $sec) = (localtime())[5,4,3,2,1,0];
my $ymdhms = sprintf("%04d%02d%02d%02d%02d%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec) ;
my $tmpfile = "$ymdhms.gb";
open (my $out, '>', $tmpfile) || die "Can't open file:$!";
my $a = 0;
my $b = $querymax-1;
while ($b<=$numtotal && $a<=$numtotal){
	my ($retmax, $retstart) = (500,0);
	#my ($retmax, $retstart) = (20,0); # debug
	print "a=$a -> b=$b\n";
	my $factory = Bio::DB::EUtilities->new(-eutil   => 'esearch',
		                               -email   => 'l.ranjard@auckland.ac.nz',
		                               -db      => 'nucleotide',
				               #-id      => \@access_num[$a..$b],
                                               -term  => join(',',@access_num[$a..$b]),
                                       	       -usehistory => 'y');
	print $factory->get_count . "\n";
	my $hist  = $factory->next_History || die 'No history data returned';
	print "History returned\n";
	# note db carries over from above
	$factory->set_parameters(-eutil   => 'efetch',
		                 -rettype => 'gb',
		                 -history => $hist);
        $retry = 0;
	RETRIEVE_SEQS:
	while (($retstart+$retmax)<=$querymax && ($a+$retstart)<=$numtotal) {
	    print "\trestart=$retstart -> retmax=$retmax\n";
	    #$factory->set_parameters(-retmax   => min($retmax,$querymax-$retstart),
	    $factory->set_parameters(-retmax   => $retmax,
		                     -retstart => $retstart);
	    eval{
		$factory->get_Response(-cb => sub {my ($data) = @_; print $out $data} );
	    };
	    if ($@) {
		die "Server error: $@.  Try again later" if $retry == 5;
		print STDERR "Server error, redo #$retry\n";
		$retry++ && redo RETRIEVE_SEQS;
	    }
	    $retstart += $retmax;
	    print "Retrieved $retstart sequences from genbank\n";
	}
	$a = $b+1;
	$b = min($b+$querymax, $numtotal);
}
close $out;

### load in the genbank file
my $seqio = Bio::SeqIO->new(-file   => $tmpfile,
                            -format => 'genbank');


#my $db = Bio::DB::Taxonomy->new(-source    => 'flatfile' ,
#                                -nodesfile => '/home/lranjard/project/nzgl01264/active/ncbi_tax/nodes.dmp',
#                                -namesfile => '/home/lranjard/project/nzgl01264/active/ncbi_tax/names.dmp');


### load NCBI taxonomy from local file
print "Retrieving taxonomy ...\n";
use Bio::LITE::Taxonomy::NCBI;
my $taxDB = Bio::LITE::Taxonomy::NCBI->new(
                                     names=>"/home/lranjard/ncbi_tax/names.dmp",
                                     nodes=>"/home/lranjard/ncbi_tax/nodes.dmp"
                                    );

my $j = 0;		# overall index of queries
my $q = 0; 		# index of query array
my $hitcount = 1;	# number of blast hits in a given query
print OUT "NoBlast" . "\t" . "unclassified;" . "\n" ; # for sequences that didn't get any Blast results, no taxonomy
while (my $seq = $seqio->next_seq){
	#print $query_nam[$q] . " " . $number_of_hits[$q] . "\n";
	while ($hitcount>$number_of_hits[$q]){
		#if ($number_of_hits[$q]==0) {
		#}
		$q++;
		$hitcount=1;
	}
	$query_nam[$q] =~ s/[;:]/_/g; # replace ";:" with "_" in the query name
	print OUT @query_nam[$q] . "_" . $hitcount . "\t";

	$hitcount++;
	### print all features of sequence object:
	#for my $feat_object ($seq->get_SeqFeatures) {          
	#   print "primary tag: ", $feat_object->primary_tag, "\n";          
	#   for my $tag ($feat_object->get_all_tags) {             
	#      print "  tag: ", $tag, "\n";             
	#      for my $value ($feat_object->get_tag_values($tag)) {                
	#	 print "    value: ", $value, "\n";             
	#      }          
	#   }       
	#}
	###
	#print $j . "\n";
	#my @organism=(); # array
	my @db_xref=(); # array
	for my $feat_object ($seq->get_SeqFeatures) {  # need to go through all features and check the one that has the tag "organism", when found put the value into the array      
	   #push @organism,$feat_object->get_tag_values("organism") if ($feat_object->has_tag("organism"));  
	   push @db_xref,$feat_object->get_tag_values("db_xref") if ($feat_object->has_tag("db_xref")); # get taxon id
	}
	#$db_xref[0] =~ s/taxon://g;
	# need to find the db_xref that has "taxon" database reference (some have BOLD database or others)
	$i = 0;
	my $taxonid = "_" ;
	foreach (@db_xref) {
		#print $db_xref[$i] ;
		if ( $db_xref[$i] =~ m/taxon:/ ) {
			$db_xref[$i] =~ s/taxon://g ;
			$taxonid = $db_xref[$i] ;
			last ;
		}
		$i++;
	}
	print $taxonid . "\n";

	#my $node = $db->get_Taxonomy_Node(-taxonid => $db_xref[0]);
	#print $node->id, " ", $node->scientific_name, " ", $node->rank, "\n";

	#my $parent = $node;
        #for (1..25) {
        #while (defined $parent->parent_id){
	#  $parent = $db->get_Taxonomy_Node(-taxonid => $parent->parent_id);
          #print "Name is ",$parent->scientific_name,"\n";
        #}
	
	if ( $taxonid =~ m/_/ ){
		print OUT "notaxonidfoundinblastresults;\n";
	} else {
		my @taxDB = ();
		@taxDB = $taxDB->get_taxonomy($taxonid);
		s/ +/_/g foreach @taxDB; # replace " " with "_" in the taxonomy
	  	#print "$_\n" for (@taxDB);
		if (!@taxDB){ # array empty
			print OUT "notaxonomyfoundinlocaldb;\n";
		} else {
	  		print OUT "$_;" for (@taxDB);
			print OUT "\n";
		}
	}

	$j++;
}

close OUT;
#unlink $tmpfile;
print "--DONE--\n";

#mothur "#classify.otu(taxonomy=uparse5.taxonomy, list=uparse5.taxonomy.list)"

