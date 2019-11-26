
#Assigment 3. María Ortiz Rodríguez
require 'bio'
require 'rest-client'

#Part 1: we take all the genes from Arabidopsis List and its sequence
def load_from_file(filename)
  file = File.open(filename, "r")
  genes = []

  file.each_line do |line|
    genes << line.delete!("\n")
  end
  file.close
  return genes
end

def seq(gene_id)
  address = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_id}"

  response = RestClient::Request.execute(
      method: :get,
      url: address)
  data = response.body
  entry = Bio::EMBL.new(data)
  seq = entry.to_biosequence

  return seq
end

#Part 2: We examine the exons and look for the CTTCTT sequence
def match_exon (seq) #object to take the target in given exon

  len_seq = seq.length()
  exon_positions = Hash.new

  # #we take all the positions inside of the seq that match with cttctt
  positions_direct_match = seq.gsub(/#{$match}/).map{Regexp.last_match.begin(0)}
  positions_reverse_match = seq.reverse_complement.gsub(/#{$match}/).map{Regexp.last_match.begin(0)}

  #Look for the features for the seq: exon id, position and strand
  seq.features.each do |feature|

    if feature.feature =~ /exon/ #first filter --> has to be exon

      position = feature.position

      exon_id = feature.qualifiers[0].value.gsub('exon_id=', '')


      if position =~ /complement/ # If that happens --> the strand is reverse -
        #coordinates for exons
        position = feature.position.match(/[0-9]+\.\.[0-9]+/).to_s.split("..")
        positions_exon = []

        position.each do |pos|
          #We transform the coordinates
          positions_exon.insert(0, len_seq - pos.to_i)
        end

        strand = '-'
        positions_reverse_match.each do |match_init|

          match_end = match_init + $match.length() - 1

          if (match_init >= positions_exon[0].to_i) && (match_init <= positions_exon[1].to_i) && (match_end >= positions_exon[0].to_i) && (match_end <= positions_exon[1].to_i)

            m_end = len_seq - match_end
            m_init = len_seq - match_init
            exon_positions[[m_end, m_init]] = [exon_id, strand]

          end

        end

        ####################################################################################################################33
      else #the exon is +

        position_exon = position.split("..")

        #target_pos_in_exon= match_exon(exon_id, positions_direct_match, len_seq, position_exon, '+')
        #
        strand = '+'
        positions_direct_match.each do |match_begin|

          match_end = match_begin + $match.length() - 1

          if (match_begin >= position_exon[0].to_i) && (match_begin <= position_exon[1].to_i) && (match_end >= position_exon[0].to_i) && (match_end <= position_exon[1].to_i)
            #The target have to be inside of the exon
            exon_positions[[match_begin, match_end]] = [exon_id, strand]

          end

        end

      end #if complement === reverse


    end #end  exon
  end #end features

    return exon_positions


end

#Part 3: We take the coordinates and create new Features
def add_features_embl(gene_id, hash_features, seq)
  # Method that iterates over the hash with the target's matched in exons
  # to add them as new features to the Bio:EMBL objects.

  exon_features = Array.new

  hash_features.each do |position, exonid_strand| ####BUENO

    feat = Bio::Feature.new("#{$match.upcase}_in_exon", "#{position[0]}..#{position[1]}")

    feat.append(Bio::Feature::Qualifier.new('nucleotide motif', "#{$match.upcase}_in_#{exonid_strand[0]}"))

    feat.append(Bio::Feature::Qualifier.new('strand', exonid_strand[1]))

    $genes_gff3.puts "#{gene_id}\t.\t#{feat.feature}\t#{position[0]}\t#{position[1]}\t.\t#{exonid_strand[1]}\t.\tID=#{exonid_strand[0]}"

    exon_features << feat
  end

  seq.features.concat(exon_features) # We add the new features created to the existing ones
end

#Part 4a: We create a GFF3 file
def chromosome_features (gene_id, seq)

  #We take the information of genes that is in primary_accession separate for ":"
  chrom_array = seq.primary_accession.split(":")


  $chromosomes_gff3.puts "#{chrom_array[2]}\t.\tgene\t#{chrom_array[3]}\t#{chrom_array[4]}\t.\t+\t.\tID=#{gene_id}"

  return [chrom_array[2], chrom_array[3], chrom_array[4]]

end

#Part 5: We create a GFF3 with chromosome coordinates
def chromosome_transform(gene, hash_features, chr)

  hash_features.each do |positions, exonid_strand|
    pos_ini_chr = chr[1].to_i + positions[0].to_i
    pos_end_chr = chr[1].to_i + positions[1].to_i

    $chromosomes_gff3.puts "#{chr[0]}\t.\tnucleotide_motif\t#{pos_ini_chr}\t#{pos_end_chr}\t.\t#{exonid_strand[1]}\t.\tID=#{exonid_strand[0]};parent=#{gene}"
  end


end

def create_file(filename)

  return File.open(filename, "w")

end


$match = "cttctt"

$genes_gff3 = create_file("genes_features.gff3")

#Part 4b: We create a file with the gene_id that doesnt have CTTCTT target
$no_match = create_file("genes_without_cttctt.txt")

$chromosomes_gff3 = create_file("chr_features.gff3")

genes = load_from_file("ArabidopsisSubNetwork_GeneList.txt")

genes.each do |gene|
  puts "Reading gene #{gene} and its sequence"
  seq_obj = seq(gene) # Create the Bio:EMBL object from each gene


  unless seq_obj.nil?
    hash_features = match_exon(seq_obj) # We get the targets inside exons of the gene


    if hash_features.empty?
      $no_match.puts gene

    else

      add_features_embl(gene, hash_features, seq_obj) # We create new features and add them to each seq_obj

      chr = chromosome_features(gene, seq_obj) # We return the chromosome number and postions

      chromosome_transform(gene, hash_features, chr) # We convert the positions to the ones that correspond in the chromosome


    end


  end

end