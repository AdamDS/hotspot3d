package TGI::Mutpro::Main::Significance;
#
#----------------------------------
# $Authors: Adam D Scott and Carolyn Lou
# $Date: 2015-09-22
# $Revision:$2015-09-22
# $URL: $
# $Doc: $ cluster significance
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Getopt::Long;
use IO::File;
use FileHandle;

sub new {
	my $class = shift;
	my $this = {};
	$this->{'hupc_file'} = undef;
	$this->{'pairwise_prefix'} = undef;
	$this->{'pairwise_suffix'} = undef;
	$this->{'clusters_prefix'} = undef;
	$this->{'clusters_suffix'} = undef;
	$this->{'AAabbrv_file'} = undef;
	$this->{'Proximity_files'} = undef;
	$this->{'simulation_number'} = 1000000;
	bless $this, $class;
	$this->process();
	return $this;
}

sub process {
	my $this = shift;
	my ( $help, $options );
	unless( @ARGV ) { die $this->help_text(); }
	$options = GetOptions (
		'hupc-file=s' => \$this->{'hupc_file'},
		'pairwise-prefix=s' => \$this->{'pairwise_prefix'},
		'ptype=s' => \$this->{'ptype'},
		'clusters-prefix=s' => \$this->{'clusters_prefix'},
		'ctype=s' => \$this->{'ctype'},
		'AAabbrv-file=s' => \$this->{'AAabbrv_file'},
		'Proximity-files=s' => \$this->{'Proximity_files'},
		'simulation-number=i' => \$this->{'simulation_number'},
		'help' => \$help,
	);
	#my $NSIMS = $this ->{'simulation_number'};

	if ( $help ) { print STDERR help_text(); exit 0; }
	unless( $options ) { die $this->help_text(); }

#if hupc file missing
	if(not defined $this->{'hupc_file'}){
		warn 'You must provide a hupc file !', "\n"; 
		die $this->help_text();
		if(not -e $this->{'hupc_file'}){
			warn "The input hupc file(".$this->{'hupc_file'}.") does not exist! ", "\n";
			die $this->help_text();
		}
	}

#if aaabbrv file missing
	unless( $this->{'AAabbrv_file'} ) {
		warn 'You must provide aaabbrv file! ', "\n";
		die $this->help_text();
	}
	unless( -e $this->{'AAabbrv_file'} ) {
		warn "The input aaabbrv file (".$this->{'AAabrv_file'}.") does not exist! ", "\n";
		die $this->help_text();
	}

#if proximity file directory missing
	unless( $this->{'Proximity_files'} ) {
		warn 'You must provide a proximity file directory ! ', "\n";
		die help_text();
	}
	unless( -d $this->{'Proximity_files'} ) {
		warn 'You must provide a valid proximity file directory ! ', "\n";
		die help_text();
	}

#processing procedure
# parse hupc file
	my ($genesOnPDB_hash_ref, $uniprot2HUGO_hash_ref, $HUGO2uniprot_hash_ref) = $this->processHUPC( $this->{'hupc_file'} );
#parse aaabbrv file
	my $AA_hash_ref = $this->processAA($this->{'AAabbrv_file'});
	$this ->getDistances($this->{'Proximity_files'}, $genesOnPDB_hash_ref, $HUGO2uniprot_hash_ref, $AA_hash_ref, $this->{'simulation_number'})
} 

# process HUPC file
sub processHUPC{
	my($this, $hupc_f) = @_;
	my $fh = new FileHandle;
	unless ($fh -> open($hupc_f)){die "Could not open hupc file ! \n"};
	my ( %genesOnPDB, %uniprot2HUGO, %HUGO2uniprot);
	while (my $line = $fh->getline){
		chomp($line);
		my ( $hugo , $uniprot , $pdb , $chain ) = split( "\t" , $line );
		$genesOnPDB{$hugo}{$pdb}{"[".$chain."]"} = 1;
		$uniprot2HUGO{$uniprot}{$hugo} = 1;
		$HUGO2uniprot{$hugo} = $uniprot;
	}
	$fh->close();
	return (\%genesOnPDB, \%uniprot2HUGO, \%HUGO2uniprot);
}

#read in amino acid information
sub processAA {
	my($this, $aaabbrv_f) = @_;
	my $fh = new FileHandle;
	unless($fh->open($aaabbrv_f)){die "Could not open AAABBRV file ~ \n"};
	my %AA;
	while ( my $line = $fh->getline ) { #file is organized with single letter, name, and three-letter abbreviation of amino acids
		chomp( $line );
		my $abbrv = (split( /\t/ , $line ))[-1];
		$AA{uc $abbrv} = 1;
	}
	$fh->close();
	return \%AA;
}

sub getDistances{
	my ($this, $proximity_d, $genesOnPDBref, $HUGO2uniprotref, $AAref, $NSIMS) = @_;
	foreach my $hugo(keys %$genesOnPDBref) {
		my $uniprot = $HUGO2uniprotref->{$hugo};
		my %distances;
		my $proxfile = $proximity_d.$uniprot."ProximityFile.csv";
		my $fh = FileHandle-> new("$proxfile", "r");
		while(my $line = <$fh>) {
			chomp $line;
			my ( $up1 , $chain1 , $res1 , $aa1 , $up2 , $chain2 , $res2 , $aa2 , $dist , $pdb , $pval ) = (split( /\t/ , $line ))[0..2,4,7..9,11,14..16];
			if ( exists $genesOnPDBref->{$hugo}{$pdb}{$chain1}
				&& exists $genesOnPDBref->{$hugo}{$pdb}{$chain2}
				#&& $chain1 eq $chain2 #specific check to mimic SpacePAC restriction
				&& exists $AAref->{$aa1}
				&& exists $AAref->{$aa2} ) {
				push @{$distances{$pdb}{$chain1}} , $dist;
			}
		}
		$fh -> close();
		
		foreach my $pdb ( keys %distances ) {
			foreach my $chain ( keys %{$distances{$pdb}} ) {
				my @distances = sort {$a <=> $b} @{$distances{$pdb}{$chain}};
				#my @distances = sort {$a <=> $b} keys %{$distances{$pdb}{$chain}};
				my $numDistances = scalar( @distances );

				#print STDOUT $tossedPairs{$pdb}{$chain}." pairs thrown out\n";
				my $fpairwise = join( "." , ( $this->{'pairwise_prefix'} , $hugo , $pdb , $chain ,$this->{'ptype'} ) );
				my $PAIRWISE = FileHandle->new( "$fpairwise" , "r" );
				if ( not defined $PAIRWISE ) {
					warn "ADSERROR: Could not open/read $fpairwise\n";
					next;
				}

				my $fclusters = join( "." , ($this->{ 'clustersPrefix'} , $hugo , $pdb , $chain ,$this->{'ctype'} ) );
				my $CLUSTERS = FileHandle->new( "$fclusters" , "r" );
				if ( not defined $CLUSTERS ) {
					warn "ADSERROR: Could not open/read $fclusters\n";
					next;
				}
				my %clusters;
				my $cluster;
				
				my $fout = $fclusters.".".$NSIMS."sims.avgDist_pvalue";
				my $OUT = FileHandle->new( "$fout" , "w" );
				if ( not defined $OUT ) {
					die "ADSERROR: Could not open/write $fout\n";
					next;
				}

				while ( my $line = <$CLUSTERS> ) {
					chomp( $line );
					my ( $clusterid , $gene , $mutation ) = (split( /\t/ , $line ))[0..2];
					$clusters{$gene}{$mutation} = $clusterid;
				}
				$CLUSTERS->close();
				
				print "Cluster\tGene\tPDB_ID\tChain\tNum_Residues\tNum_Pairs\tAvg_Distance\tEstimated_PValue\n";
				$OUT->print( "Cluster\tGene\tPDB_ID\tChain\tNum_Residues\tNum_Pairs\tAvg_Distance\tEstimated_PValue\n" );

				my %mass2pvalues;
				my %mapping;
				my %withinClusterSumDistance;
				while ( my $line = <$PAIRWISE> ) {
					chomp( $line );
					my @line = split( /\t/ , $line );
					my ( $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 , $info ) = @line[0,4,5,6,9,13,14,15,19];
					if ( exists $clusters{$gene1}{$mu1}
						&& exists $clusters{$gene2}{$mu2}
						&& $clusters{$gene1}{$mu1} == $clusters{$gene2}{$mu2} ) {
						my @infos = split( /\|/ , $info );
						foreach my $dinfo ( @infos ) {
							my ( $dist , $pdb , $pval ) = split( /\s/ , $dinfo );
							if ( exists $genesOnPDBref->{$gene1}{$pdb}{$chain1}
							&& exists $genesOnPDBref->{$gene2}{$pdb}{$chain2}
							&& $chain1 eq $chain2 ) { #specific check to mimic SpacePAC restriction
								$cluster = $clusters{$gene1}{$mu1};
								$mapping{$cluster}{$res1} = 1;
								$mapping{$cluster}{$res2} = 1;
								$withinClusterSumDistance{$cluster} += $dist;
								$mass2pvalues{$cluster} = $pval;
							} #if gene,pdb,chain1 exists
						} #foreach dinfo
					} #if genes/mus exists in clusters
				} #while pairwise file
				$PAIRWISE->close();
	
				my %massDistribution;
				my %withinClusterAvgDistance;
				my ( $mass , $numPairs );
				foreach $cluster ( keys %mapping ) {
					$mass = scalar keys %{$mapping{$cluster}};
					$numPairs = ( $mass ) * ( $mass - 1 ) / 2;
					if ( $mass > 2 ) {
						$massDistribution{$mass}{$cluster} = 1;
						$withinClusterAvgDistance{$cluster} = $withinClusterSumDistance{$cluster} / $numPairs;
						delete $mass2pvalues{$cluster};
					} elsif ( $mass == 2 ) {
#print STDOUT "Cluster ".$cluster." has ".$mass." residues, so it's p-value is ".$mass2pvalues{$cluster}."\n";
						$OUT->print( join( "\t" , ( $cluster , $hugo , $pdb , $chain , $mass , $numPairs , $withinClusterSumDistance{$cluster} , $mass2pvalues{$cluster} ) )."\n" );
					} elsif ( $mass == 1 ) {
						#print STDOUT "Cluster ".$cluster." has ".$mass." residue, so it has no p-value\n";
						$OUT->print( join( "\t" , ( $cluster , $hugo , $pdb , $chain , $mass , $numPairs , 0 , "NA" ) )."\n" );
					} #if mass block
				} #foreach cluster

				#my %permutationTestPValues;
				#print STDOUT $hugo."\t".$pdb."\t".$chain."\n";
				foreach my $numResidues ( keys %massDistribution ) {
					$numPairs = ( $numResidues ) * ( $numResidues - 1 ) / 2;
					my %below;
					#print STDOUT "Simulation\tSum_Distances\tAvg_Distance\n";
					for ( my $simulation = 0; $simulation < $NSIMS ; $simulation++ ) {
						my $sumDists = 0;
						for ( my $i = 0; $i < $numPairs ; $i++ ) {
							my $randIndex = int( rand( $numDistances ) );
							$sumDists += $distances[$randIndex];
						} #foreach random pick
						my $avgDist = $sumDists / $numPairs;
						#print STDOUT $simulation."\t".$sumDists."\t".$avgDist."\n";
						foreach $cluster ( keys %{$massDistribution{$numResidues}} ) {
							if ( $avgDist < $withinClusterAvgDistance{$cluster} ) {
							$below{$cluster}++;
							}
						} #foreach cluster
					} #foreach simulation
					foreach $cluster ( keys %{$massDistribution{$numResidues}} ) {
						my $permutationTestPValue = $below{$cluster} / $NSIMS;
						#print STDOUT "Estimated p-value: ".$permutationTestPValue."\n";
						$OUT->print( join( "\t" , ( $cluster , $hugo , $pdb , $chain , $numResidues , $numPairs , $withinClusterAvgDistance{$cluster} , $permutationTestPValue ) )."\n" );
					} #foreach cluster
				} #foreach numResidues
				$OUT->close();
			} #foreach chain
		} #foreach pdb
	} #foreach gene/protein from hupc
	return 1;
}

sub help_text{
	my $this = shift;
	return <<HELP

Usage: hotspot3d significance [options]

	--maf-file				Input MAF file
	--data-dir				HotSpot3D preprocessing results directory

	--drugport-file			Drugport database parsing results file ( optional )
	--output-prefix			Prefix of output files, default: 3D_Proximity 
	--skip-silent			Skip silent mutations, default: no
	--missense-only 		Missense mutation only, default: no
	--p-value				p_value cutoff(<=), default: 0.05
	--3d-dis				3D distance cutoff (<=), default: 10
	--linear-dis			Linear distance cutoff (>=), default: 20 
	--transcript-id-header	MAF file column header for transcript id's, default: transcript_name
	--amino-acid-header		MAF file column header for amino acid changes, default: amino_acid_change 

	--help					This message

HELP

}

1;

