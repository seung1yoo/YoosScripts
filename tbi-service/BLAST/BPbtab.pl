#!/usr/local/bin/perl


=head1 NAME
    
    BPbtab

    script parses WU-BLAST or NCBI-BLAST output files into BTAB format where each HSP is reported as a single line with tab-delimited fields.

    The original BTAB program is described here:

    Dubnick M. (1992) Btab--a Blast output parser. Comput Appl Biosci 8(6):601-2

    
    A Perl version which emulates the functionality of btab is provided here by BPbtab.  BPbtab relies exclusively on the BioPerl 1.4 Bio::SearchIO functionality.


    
=cut


=head1 USAGE

    Standard input is parsed and written to standard output.

    BPbtab blast.output > blast.output.btab

=cut


use strict;
use lib "/BiOfs/BioPeople/sgpark/BioTools/BioPerl/lib/perl5/";
use Bio::SearchIO;

my $in = new Bio::SearchIO(-format => 'blast', 
			   -file   => $ARGV[0] );

print "#q_name\tq_length\tq_start\tq_end\tq_strand\td_name\td_start\td_end\td_strand\tidentify\tbit_score\te_value\tq_cov\td_cov\tdescription\n";

# parse each blast record:
while( my $result = $in->next_result ) {

    # parse each hit per record.
    while( my $hit = $result->next_hit ) {
	
	# a hit consists of one or more HSPs
	while ( my $hsp = $hit->next_hsp ) {
	    my @x;
	    $x[0] = $result->query_name();
	    $x[1] = $result->query_length();
	    $x[2] = $hsp->start('query');
	    $x[3] = $hsp->end('query');

	    my $queryStrand = $hsp->strand('query');
            my $strandDescript = "null";
            if ($queryStrand == 1) {
                $strandDescript = "Plus";
            } elsif ($queryStrand == -1) {
                $strandDescript = "Minus";
            }

	    $x[4] = $strandDescript;
	    $x[5] = $hit->name();
	    $x[6] = $hsp->start('hit');
	    $x[7] = $hsp->end('hit');
	    my $hitStrand = $hsp->strand('hit');
	    
	    $strandDescript = "null";
	    if ($hitStrand == 1) {
		$strandDescript = "Plus";
	    } elsif ($hitStrand == -1) {
		$strandDescript = "Minus";
	    }
	    $x[8] = $strandDescript;
#	    if ($x[8] =~ /Plus/) { next; }
	    $x[9] = sprintf ("%.1f", $hsp->percent_identity());   
#	    if ( $x[9] < 95 ) { next; }
#	    if ( (($x[3]-$x[2]+1)/$x[1]) < 0.90 && (($x[7]-$x[6]+1)/$hit->length()) < 0.90 ) { next; }
#	    if ( $x[0] =~ /^$x[5]$/) { next; }
#	my @name = split(/[_-]/,$x[0]);
#	unless ($name[1] =~ /^$x[5]$/) { next; }

	    $x[10] = $hsp->score();
	    $x[11] = $hsp->evalue();
	    $x[12] = sprintf ("%.1f", ($x[3]-$x[2]+1)/$x[1]*100);
	    $x[13] = sprintf ("%.1f", (($x[7]-$x[6]+1)/$hit->length()*100));
	    if ($x[12] < 50 and $x[13] < 50) { next; }
	    $x[14] = $hit->description;   
#	    $x[11] = $hsp->query_string();
#	    $x[13] = $hsp->hit_string();
	    my $outline = join ("\t", @x);
	    print "$outline\n";
	}
    }
}


=head1 DESCRIPTION

                                   Fields


The parsed BLAST output is presented as single line HSP descriptions, with tab-delimited fields in the following order:


[0]      Query Sequence Name

[2]      Query Sequence Length

[3]      Search Method  --  Blast family application name

[4]      Database Name

[5]      Subject Sequence Name  --  Database entry name

[6],[7]   Query Left End, Query Right End  --  The endpoints of the part
               of the query sequence which Blast aligns with the subject sequence.

[8],[9]  Subject Left End, Subject Right End  --  The endpoints of the part 
               of the subject sequence which Blast aligns with the query sequence.

[10]     Percent Identity  --  The fraction of residues which are absolute
            matches between the query and subject sequence, expressed in
            percent.  

[11]     Percent Similarity  --  The fraction of residues which are exact or 
            similar matches between the query and subject sequence, expressed in
            percent. 

[12]    HSP score    

[13]    Bits score


[15]    Description  --  A freeform text field which contains the biological
        description field from the database for the subject sequence.  If
        this text occupies more than one line in the Blast output file, the
        NewLines are replaced by spaces.  Commas may occur in this field even
        if they are the field separator character, because this is the last
        field in the record.

[16]    Query Frame (1, 2, 3, -1, -2, -3)

[17]    Query Strand  --  Plus, Minus or null

[18]    DB sequence length

[19]    Expect -- expected value

[20]    P-Value  --  Poisson ratio


** Note ** Intervening field positions which are not described are not currently supported.  These remain to support compatibility with other existing btab implementations.

=cut
