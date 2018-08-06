#!/usr/bin/perl

use strict;
use warnings;
use autodie;

use Carp;
use Getopt::Long;
#use List::Util      qw( first  );
#use List::MoreUtils qw( notall );
#use File::Basename;
use Bio::SeqIO;

my $input_sequence;
my $input_gff;
my $new_id;
my $output_sequence;
my $output_gff;
my $show_usage;
my $show_help;
my $show_man;

my $format = q{fasta};

my $usage = qq{Usage:\n\n$0 [ -sequencefile ] SEQUENCE_FILE [ [ -gff ] GFF_FILE ] [ [ -newid ] NEW_SEQ_ID ] [options]};

my $help = qq{$usage

Options:

	-sequencefile	STRING
	-format		STRING
	-gff		STRING
	-newid		STRING
	-newsequence	STRING
	-newgff		STRING
};

my $man = qq{$usage

Options:

	-sequencefile	STRING	name of input sequence file

	-format		STRING	format of input sequence file (and output
				sequence file); must be supported by
				Bio::SeqIO. Default is '$format' .

	-gff		STRING	name of input GFF file

	-newid		STRING	identifier ('accession') of the single
				sequence which will be output (to both the
				output sequence file and GFF file); default
				is the same as the ID of the first sequence
				encountered

	-newsequence	STRING	name of output sequence file containing a
				single sequence record; if omitted, then the
                                script will still run, and commentate on
                                sequence lengths, offsets etc;

	-newgff		STRING	name of output GFF file containing annotations
				for a single sequence; can also be omitted


Purpose:

The background to this script is that if you want to use genomic sequence
files and GFF (or other format) annotation files with Artemis or ACT, things
can get quite messy if your sequence data is split across more than one file.
That is, if you have genome data which comes in the form of several or many
contig/scaffold sequences, you really need to concatenate this into a single
sequence (there are facilities to do this within Artemis/ACT). This single
sequence of course has continuous coordinates (i.e. base-pair numbering
beginning at 1 and incrementing across all the previous inter-contig
'boundaries', which are no longer there). In effect, all but the first
sequence has its coordinates renumbered.

If you then try to import annotations (from GFF etc), all coordinates of
features will therefore be wrong, apart from the features of the first
sequence.

Note that simply loading a GenBank-format file (which might include feature
annotations as well as the sequences) containing multiple sequences, i.e.
what is provided by for example NCBI GenBank and RefSeq, does not
seem to work either because Artemis appears to assume that this pertains to
only a single sequence, and stops after it has read the data of the first
sequence; at least, that is all that is displayed, and it provides a list
of only 1 'entry', so there is no way of getting at sequences 2, 3, .. etc.

If the sequence and annotation data is in one file per contig/scaffold, then
reading in sequence data one file at a time is also not the way. It appears
that once the first sequence is in, the full-range coordinate system is
defined, and subsequent reading of entries appears to have to fit this system.

For example, if you read in the first contig sequence, which is 6789 bp long,
then you have a coordinate system of 1 .. 6789. If you then try to read in
a second entry (from the first entry's File -> Read An Entry ... dialogue;
not the principal Artemis window's File --> Open ... option, which does not
provide the solution either) - which is longer than 6789 bp - then an error
is produced, about 'out of range location'.

Therefore, the simplest solution seems to be - to not just concatenate the
sequence file, but to concatenate the GFF entries into a single sequential
set of annotations - with the begin .. end coordinates all adjusted
appropriately.

This script therefore concatenates a sequence file, and optionally
concatenates a GFF file - while changing the sequence IDs and begin, end
coordinates in each GFF line, in accordance with the new, single-sequence
FASTA file.

Output of a new FASTA and/or GFF file is optional; if specified or not, the
script still reports which changes would be (or are being) made.

};

croak $usage if !GetOptions(

    'format=s'       => \$format,
    'sequencefile=s' => \$input_sequence,
    'gff=s'          => \$input_gff,
    'newid=s'        => \$new_id,
    'newsequence=s'  => \$output_sequence,
    'newgff=s'       => \$output_gff,
    'usage'          => \$show_usage,
    'help'           => \$show_help,
    'man'            => \$show_man,
   
);

if ($show_usage) { print $usage; exit 0; }
if ($show_help)  { print $help;  exit 0; }
if ($show_man)   { print $man;   exit 0; }

$input_sequence ||= shift;
$input_gff      ||= shift;
$new_id         ||= shift;

#croak qq{input sequence and annotation (GFF) files are mandatory}
#    if !$input_sequence || !$input_gff;

croak qq{input sequence is mandatory} if !$input_sequence;

my $seqio = Bio::SeqIO->new(
    '-file'   => $input_sequence,
    '-format' => $format,
);

my $coord_offset = 0;

my %offset_of_seq;
my %length_of_seq;

my $concatenated_sequence = q{};
my $concatenated_description = q{concatenation of sequences:}; # summarises the IDs, not the actual descriptions

while ( my $seq_obj = $seqio->next_seq() ) {

    my $seq_id     = $seq_obj->display_id();
    my $seq_length = $seq_obj->length();
    $concatenated_sequence .= $seq_obj->seq();
    $concatenated_description .= qq{ $seq_id ($seq_length bp)};

    printf qq{sequence %-30s\t%8d bp\tcoordinates will be offset by +%8d bp\n},
        $seq_id, $seq_length, $coord_offset;
    $offset_of_seq{$seq_id} = $coord_offset;
    $length_of_seq{$seq_id} = $seq_length;
    $coord_offset += $seq_length;

    $new_id ||= $seq_id;
}

print qq{total sequence length: $coord_offset bp\n};

if ($output_sequence) {
    print qq{creating new sequence file '$output_sequence', containing 1 sequence '$new_id' of $coord_offset bp\n};
    my $seq_obj = Bio::Seq->new( '-display_id' => $new_id,
                                 '-seq' => $concatenated_sequence,
                                 '-desc' => $concatenated_description);
    my $seqio = Bio::SeqIO->new(
        '-file'   => ">$output_sequence",
        '-format' => $format,
    );
    $seqio->write_seq($seq_obj);
};

exit if !$input_gff;

# really no need to use an API to read the GFF file

open my  $fh, '<', $input_gff;

open my $ofh, '>', $output_gff if $output_gff;

LINE:
while (defined (my $line = <$fh>)) {

    my $new_line = $line;

    ($line =~ m{ \A [#][#] sequence [-] region \s+ (\S+) \s+ (\d+) \s+ (\d+) }xms) && do {
        my ($seq_id, $old_begin, $old_id) = ($1, $2, $3);
        croak qq{unknown sequence '$seq_id'} if !defined $offset_of_seq{$seq_id};
        # only make a change if the offset is non-zero, or if the ID needs
        # changing (in most cases, one or both of these will be true)
        if ($offset_of_seq{$seq_id} || ($seq_id ne $new_id)) {
            my $new_begin = 1 + $offset_of_seq{$seq_id};
            my $new_end   = $offset_of_seq{$seq_id} + $length_of_seq{$seq_id};
            $new_line  = qq{##sequence-region $new_id $new_begin $new_end\n};
            print qq{changing:\n${line}to:\n$new_line\n};
        }
        print $ofh $new_line if defined $ofh;
        next LINE;
    };

    ($line =~ m{ \A [#] }xms) && do {
        print $ofh $line if defined $ofh;
        next LINE;
    };

    my ($seq_id, $source, $record_type, $begin, $end, $remainder)
       = split m{ \t }xms, $line, 6;

    croak qq{unknown sequence '$seq_id'} if !defined $offset_of_seq{$seq_id};

    if ($offset_of_seq{$seq_id} || ($seq_id ne $new_id)) { # i.e. non-zero offset or new ID
        my $new_begin = $begin + $offset_of_seq{$seq_id};
        my $new_end   = $end   + $offset_of_seq{$seq_id};
        $new_line = join(qq{\t}, $new_id, $source, $record_type,
                            $new_begin, $new_end, $remainder);
        print qq{changing:\n${line}to:\n$new_line\n};
    }

    print $ofh $new_line if defined $ofh;
}

close $ofh if defined $ofh;
close $fh;
