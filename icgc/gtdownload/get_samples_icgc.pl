use strict;
use XML::DOM;
use Data::Dumper;
use Getopt::Long;
use XML::LibXML;
use Cwd;
use URI::Split qw( uri_split uri_join );
use String::Util qw(trim);
use Storable;

# DESCRIPTION
# A tool for identifying samples for download from centers  
# TODO:
# * need to use perl package for downloads, not calls out to system

#############
# VARIABLES #
#############

my %hash_downloaded_files = ();
my $gnos_url = "https://gtrepo-ebi.annailabs.com";
# download data by default 
my $download_files = 1;
my $cgquery='/usr/bin/cgquery';
my $gtdownload='/usr/bin/gtupload';
my $no_download_results=0;
my $dcc_project_list = 'dcc_projects.txt';
# dir for ourtputing the files
my $output_dir = cwd();
# keep track of samples downloaded from each data centre
my $report_name = "report_samples.txt";
my $alicot_ids = {};
my $cntrl_ids = {};
# or you can define them as 
# my %alicot_ids = ();
# and fill it as:
# $alicot_ids{$a} = $b;
# otherwise
# my $alicot_ids = {};
# then
# $alicot_ids->{$a} = $b;

my $center = get_center_info($gnos_url);
my $hashfile = "files.".$center.".hash";

if (scalar(@ARGV) < 1 || scalar(@ARGV) > 2) {
  #print "USAGE: 'perl $0 --gnos-url <URL> [--working-dir <working_dir>] [--report <download_report.txt>] [--no-download]'\n";
  print "USAGE: 'perl $0 --gnos-url <URL> [--dcc_projects <dcc_projects_list.txt>] [--no-download]'\n";  
  print "\t--gnos-url a URL for a GNOS server, e.g. https://gtrepo-ebi.annailabs.com\n";
  print "\t--dcc-projects provide a new list of dcc project codes\n";
  print "\t--output-dir for downloading the files (cwd by default)\n";
  print "\t--report the report file name\n";
  print "\t--no-download do not download files (false by default)\n";
  exit;
}

GetOptions("gnos-url=s" => \$gnos_url, "dcc-projects:s" => \$dcc_project_list, "no-download" => \$no_download_results, "report_name=s" => \$report_name, "output_dir=s" => \$output_dir );

if ($no_download_results) { $download_files = 0; }

##############
# MAIN STEPS #
##############

# get the list of dcc projects 
open(my $fh, '<', $dcc_project_list)
    or die("Can't open input file \"$dcc_project_list\": $!\n");
my %dcc_codes;
     while (my $line=<$fh>) {   
     chomp;
   my ($code,$name) = split /\t/, $line;  
          $dcc_codes{trim($code)} = trim($name);   
    }

# check dcc_codes hash
#print Dumper \%dcc_codes;

# output report file
open R, ">$report_name" or die;

# READ CLUSTER INFO AND RUNNING SAMPLES
#my ($cluster_info, $running_samples) = read_cluster_info($cluster_json);
#print Dumper($cluster_info);
#print Dumper($running_samples);

# READ INFO FROM GNOS
#my $sample_info = read_sample_info();
#print Dumper($sample_info);

# setup directories donwload samples ...

if (-e "$output_dir/$hashfile")
{
 print "FOUND HASH\n";
 print "LOADING HASH\n";   
 %hash_downloaded_files = %{retrieve($hashfile)};
 print "size of hash:  " . keys( %hash_downloaded_files ) . ".\n";
 #print Dumper \%hash_downloaded_files;
} 

# pass a hash to the subroutine 
#get_samples(\%hash_downloaded_files);
get_samples();


close R;

open CNTRL, ">report_controls.txt" or die;

my $total_file_size_cntrl = 0;

print "DOWNLOADING MATCHED CONTROLS ...\n";

# DOWNLOAD MATCHED CONTROLS 
foreach my $cntrl(keys %{$cntrl_ids}) {
    print CNTRL "CNTRL for SAMPLE: $cntrl\n";
    my $num = $alicot_ids->{$cntrl_ids->{$cntrl}}{sample};
    #my $origin_tumour = $alicot_ids->{$cntrl_ids->{$cntrl}}{analysis};
    print CNTRL "$num\n";
    #print CNTRL "$origin_tumour\n";
    my $size_gb = get_cntrl_sample($num);
    $total_file_size_cntrl = $total_file_size_cntrl + $size_gb;
}

print CNTRL "BAM total file size: $total_file_size_cntrl GB\n";

close CNTRL;

# print Dumper \%hash_downloaded_files;

# save hash of downloaded files 
store(\%hash_downloaded_files, $hashfile); 



###############
# SUBROUTINES #
###############

sub get_samples {
  #my $hash= shift;
  #print Dumper($hash);
  open OUT, ">xml_parse.log" or die;
  ##my $d = {};  
  # PARSE XML
  my $parser = new XML::DOM::Parser;
  my $auth = get_center_info($gnos_url);
  ##my $cmd = "mkdir -p xml/$auth; $cgquery -s $gnos_url -o xml/$auth/data.xml 'study=*&state=live'"; print OUT "$cmd\n"; system($cmd); 
  my $doc = $parser->parsefile("xml/$auth/data.xml");  
  # print OUT all HREF attributes of all CODEBASE elements
  my $nodes = $doc->getElementsByTagName ("Result");
  my $n = $nodes->getLength;
  print "$n\n";

  # DEBUG
  #my $n = 1;

  print OUT "\n";
  
  my $cnt_samples = 0;
  my $total_file_size=0;

  for (my $i = 0; $i < $n; $i++)
  {
      my $node = $nodes->item ($i);
      my $aurl = getVal($node, "analysis_full_uri");
      if($aurl =~ /^(.*)\/([^\/]+)$/) {
      $aurl = $1."/".lc($2);
      } else {
        print OUT "SKIPPING!\n";
        next;
      }
      print OUT "ANALYSIS DETAIL URL: $aurl\n";
      download($aurl, "xml/$auth/data_$i.xml");
      my $adoc = $parser->parsefile ("xml/$auth/data_$i.xml");
      my $adoc2 = XML::LibXML->new->parse_file("xml/$auth/data_$i.xml");
      my $analysisId = getVal($adoc, 'analysis_id');
      my $analysisDataURI = getVal($adoc, 'analysis_data_uri');
      my $center_name = getVal($adoc, 'center_name');
      #my $submitterAliquotId = getCustomVal($adoc2, 'submitter_aliquot_id,submitter_sample_id');
      #my $aliquotId = getCustomVal($adoc2, 'aliquot_id');
      my $aliquotId = getVal($adoc, 'aliquot_id');
      $alicot_ids->{$aliquotId}{sample} = $i;
      #my $submitterParticipantId = getCustomVal($adoc2, 'submitter_participant_id,submitter_donor_id');
      my $participantId = getCustomVal($adoc2, 'participant_id,submitter_donor_id');
      my $dcc_code = getCustomVal($adoc2,'dcc_project_code');
      my ($dcc_disease, $country) = split /-/, $dcc_code; 
      my $dcc_specimen_type = getCustomVal($adoc2,'dcc_specimen_type');
      my $disease_abbr = getCustomVal($adoc2,'disease_abbr');
      my $submitterSampleId = getCustomVal($adoc2, 'submitter_sample_id');
      # if donor_id defined then dealing with newer XML
      if (defined(getCustomVal($adoc2, 'submitter_donor_id')) && getCustomVal($adoc2, 'submitter_donor_id') ne '') {
        $submitterSampleId = getCustomVal($adoc2, 'submitter_specimen_id');
      }
      my $sampleId = getCustomVal($adoc2, 'sample_id,submitter_specimen_id');
      my $use_control = getCustomVal($adoc2, "use_cntl");
      $alicot_ids->{$use_control}{analysis} = $analysisId;
      if ($use_control eq '' || $use_control eq 'N/A' || $use_control eq 'NA') {$use_control = 'NA'};
      my $alignment = getVal($adoc, "refassem_short_name");
      #my $total_lanes = getCustomVal($adoc2, "total_lanes");
      my $sample_uuid = getXPathAttr($adoc2, "refname", "//ANALYSIS_SET/ANALYSIS/TARGETS/TARGET/\@refname");
      my $unaligned_reads = getVal($adoc, 'alignment_includes_unaligned_reads');
      #print OUT "ANALYSIS: $analysisDataURI \n";      
      #print OUT "ANALYSISID: $analysisId\n";
      #print OUT "CENTER NAME: $center_name\n";
      #print OUT "PARTICIPANT ID: $participantId\n";
      #print OUT "SAMPLE ID: $sampleId\n";
      #print OUT "ALIQUOTID: $aliquotId\n";
      #print OUT "SUBMITTER PARTICIPANT ID: $submitterParticipantId\n";
      #print OUT "SUBMITTER SAMPLE ID: $submitterSampleId\n";
      #print OUT "SUBMITTER ALIQUOTID: $submitterAliquotId\n";
      #print OUT "ALIGNMENT: $alignment\n";
      #print OUT "Control: $use_control\n"; 
      #print OUT "DCC CODE: $dcc_code\n"; 
      #print OUT "DCC DESEASE: $dcc_disease\n"; 
      #print OUT "DESEASE: $disease_abbr\n"; 
      #print OUT "DCC SPECIMEN TYPE: $dcc_specimen_type\n"; 
      #print OUT "UNALIGNED READS INCLUDED: $unaligned_reads\n";
      my $libName = getVal($adoc, 'LIBRARY_NAME');
      my $libStrategy = getVal($adoc, 'LIBRARY_STRATEGY');
      my $libSource = getVal($adoc, 'LIBRARY_SOURCE');
      my $intrument = getVal($adoc, 'INSTRUMENT_MODEL');
      my $aligner = getVal($adoc, 'NOTES');
      #print OUT "ALIGNER: $aligner\n"; 
      #print OUT "LibName: $libName LibStrategy: $libStrategy LibSource: $libSource Intrument: $intrument\n";
      # get the BAM file
      my $files = readFiles($adoc);
      #foreach my $file(keys %{$files}) {
      #   print OUT "LOCAL FILE PATH TO BAM: $analysisId/$file\n";
      #  my $size_gb = $files->{$file}{size}/(1.024e+09);
      #	 print OUT "BAM file size: $size_gb GB\n";
      #}
      # test if disease is included in DCC disease codes
      my $known_dcc_code = 0;
      if (exists $dcc_codes{$dcc_code}) {$known_dcc_code = 1}
      my $substr1 = 'tumour';
      if (!exists $hash_downloaded_files{$analysisDataURI}){
	  if($download_files && $known_dcc_code && $alignment ne 'unaligned' && index(lc($dcc_specimen_type), $substr1) != -1) {
	      #print "DOWNLOAD\n";
	      $cnt_samples++;
	      $cntrl_ids->{$i} = $use_control;
	      print R "Sample: $i\n";
	      print R "ANALYSIS: $analysisDataURI \n";
	      #print R "ALICOTEid: $alicot_ids->{$aliquotId}{analysis}\n";
	      #print R "ALICOTEid: $aliquotId\n";
	      print R "DCC CODE: $dcc_code\n";
	      print R "CENTER NAME: $center_name\n";
	      print R "DCC DESEASE: $dcc_disease\n"; 
	      print R "DCC SPECIMEN TYPE: $dcc_specimen_type\n";
	      print R "Control: $use_control\n";
	      print R "SUBMITTER SAMPLE ID: $submitterSampleId\n";
	      print R "ALIGNMENT: $alignment\n";
	      print R "ALIGNER: $aligner\n";
	      print R "LibName: $libName LibStrategy: $libStrategy LibSource: $libSource Intrument: $intrument\n";
	      foreach my $file(keys %{$files}) {
		  print R "PATH TO BAM TO DOWNLOAD: $analysisId/$file\n";
		  my $size_gb = $files->{$file}{size}/(1.024e+09);
		  print R "BAM file size: $size_gb GB\n";
		  $total_file_size = $total_file_size + $size_gb;
	          print R "FILE: $file\n";
		  print R "gtdownload -c pubkey -v -d $analysisDataURI\n";
		  print R "OUTPUT DIR: -p $output_dir/$dcc_code/$analysisId/tumour/\n";
		  my $cmd2 = "$gtdownload -c keyfile -v -d $analysisDataURI -p $output_dir/$dcc_code/$analysisId/tumour/";
		  system($cmd2);
	      }
	      $hash_downloaded_files{$analysisDataURI} = 'tumour'; 
	      #print R "gtdownload -c pubkey.txt -v -d $analysisDataURI/$file\n";
	      #-d /.mounts/labs/ferr_lab/scratch/data/TCGA/coad/coad.xml -p /.mounts/labs/ferr_lab/scratch/data/TCGA/coad/
	      #my $cmd2 = "$gtdownload -c pubkey.txt -v -d $analysisDataURI\n";
	      #system($cmd2);
	  }
      #} else {
      #	  print "FILE EXISTS WILL NOT DOWNLOAD\n";
      }

  }

  print OUT "Total number of samples downloaded: $cnt_samples\n";
  print OUT "Total file size downloaded: $total_file_size\n";  
  # Print doc file
  #$doc->printToFile ("out.xml");
  
  # Print to string
  #print OUT $doc->toString;
  
  # Avoid memory leaks - cleanup circular references for garbage collection
  #$doc->dispose;
  close OUT;
  #return($d);
}


sub get_cntrl_sample {
    my ($num) = @_;
    my $parser = new XML::DOM::Parser;
    my $auth = get_center_info($gnos_url);
    my $adoc = $parser->parsefile ("xml/$auth/data_$num.xml");
    my $adoc2 = XML::LibXML->new->parse_file("xml/$auth/data_$num.xml");
    my $analysisId = getVal($adoc, 'analysis_id');
    my $analysisDataURI = getVal($adoc, 'analysis_data_uri');
    my $alignment = getVal($adoc, "refassem_short_name");
    my $dcc_code = getCustomVal($adoc2,'dcc_project_code');
    my $aliquotId = getVal($adoc, 'aliquot_id');
    my $dcc_specimen_type = getCustomVal($adoc2,'dcc_specimen_type');
    my $libName = getVal($adoc, 'LIBRARY_NAME');
    my $libStrategy = getVal($adoc, 'LIBRARY_STRATEGY');
    my $libSource = getVal($adoc, 'LIBRARY_SOURCE');
    my $intrument = getVal($adoc, 'INSTRUMENT_MODEL');
    my $files = readFiles($adoc);
    my $size_gb = 0;
    my $substr1 = 'normal';
    my $tumor_origin = $alicot_ids->{$aliquotId}{analysis};
    if (!exists $hash_downloaded_files{$analysisDataURI}){
	if($alignment ne 'unaligned' && index(lc($dcc_specimen_type), $substr1) != -1) {
	    foreach my $file(keys %{$files}) {
	      print CNTRL "PATH TO BAM TO DOWNLOAD: $analysisId/$file\n";
	      $size_gb = $files->{$file}{size}/(1.024e+09);
	      print CNTRL "BAM file size: $size_gb GB\n";
	      print CNTRL "FILE: $file\n";
	      print CNTRL "TUMOR ORIGIN: $tumor_origin\n";
	      print CNTRL "gtdownload -c pubkey -v -d $analysisDataURI\n";
	      print CNTRL "OUTPUT DIR: -p $output_dir/$dcc_code/$tumor_origin/normal/\n";
	      my $cmd3 = "$gtdownload -c keyfile -v -d $analysisDataURI -p $output_dir/$dcc_code/$tumor_origin/normal/";
	      system($cmd3);
	  }
	$hash_downloaded_files{$analysisDataURI} = 'cntrl';
     }
    }
    
return($size_gb)

}



sub readFiles {
  my ($d) = @_;
  my $ret = {};
  my $nodes = $d->getElementsByTagName ("file");
  my $n = $nodes->getLength;
  for (my $i = 0; $i < $n; $i++)
  {
    my $node = $nodes->item ($i);
my $currFile = getVal($node, 'filename');
my $size = getVal($node, 'filesize');
my $check = getVal($node, 'checksum');
            $ret->{$currFile}{size} = $size;
            $ret->{$currFile}{checksum} = $check;
  }
  return($ret);
}



sub get_center_info {
    ## get the data center name from the http address
    my ($url) = @_;
    my ($scheme, $auth, $path, $query, $frag)  = uri_split( $url );
    $auth =~ s/gtrepo-//s; 
    $auth =~ s/.annailabs.com//s;
    return($auth);
}


sub getCustomVal {
  my ($dom2, $keys) = @_;
  my @keys_arr = split /,/, $keys;
  for my $node ($dom2->findnodes('//ANALYSIS_ATTRIBUTES/ANALYSIS_ATTRIBUTE')) {
    my $i=0;
    for my $currKey ($node->findnodes('//TAG/text()')) {
      $i++;
      my $keyStr = $currKey->toString();
      foreach my $key (@keys_arr) {
        if ($keyStr eq $key) {
          my $j=0;
          for my $currVal ($node->findnodes('//VALUE/text()')) {
            $j++;
            if ($j==$i) {
              return($currVal->toString());
            }
          }
        }
      }
    }
  }
  return("");
}

sub getXPathAttr {
  my ($dom, $key, $xpath) = @_;
  #print "HERE $dom $key $xpath\n";
  for my $node ($dom->findnodes($xpath)) {
    #print "NODE: ".$node->getValue()."\n";
    return($node->getValue());
  }
  return "";
}

sub getVal {
  my ($node, $key) = @_;
  #print "NODE: $node KEY: $key\n";
  if ($node != undef) {
    if (defined($node->getElementsByTagName($key))) {
      if (defined($node->getElementsByTagName($key)->item(0))) {
        if (defined($node->getElementsByTagName($key)->item(0)->getFirstChild)) {
          if (defined($node->getElementsByTagName($key)->item(0)->getFirstChild->getNodeValue)) {
           return($node->getElementsByTagName($key)->item(0)->getFirstChild->getNodeValue);
          }
        }
      }
    }
  }
  return(undef);
}

sub download {
  my ($url, $out) = @_;

  my $r = system("wget -q -O $out $url");
  if ($r) {
$ENV{PERL_LWP_SSL_VERIFY_HOSTNAME}=0;
    $r = system("lwp-download $url $out");
    if ($r) {
print "ERROR DOWNLOADING: $url\n";
exit(1);
    }
  }
}
