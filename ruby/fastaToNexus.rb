#! /usr/bin/ruby

# Converts a FASTA file into a NEXUS file

require_relative 'modules/FastaReader'

if ARGV.length != 2 then
    puts "[USAGE]: fastaToNexus.rb <infile.fas> <outfile.nex>"
    Process.exit
end

infile = ARGV[0]
outfile = ARGV[1]
sequences = {}
duplicates = Hash.new(1)
length = 0
irregular_in_head = false
irregulars =['-', '.', ' ', '/', '|']

fasta = FastaReader.new(infile)

fasta.parse_records do |header, sequence|
    if sequences.has_key? header then
        duplicates[header] += 1
    end
    if length == 0 then
        length = sequence.length
    else
        if sequence.length != length then
            puts "[ERROR]: Mismatch in sequence lengths: #{length} VS #{sequence.length}"
            Process.exit
        end
    end
    if irregular_in_head == false then
        if irregulars.any? { |irr| header.include? irr } then
            irregular_in_head = true
        end
    end
    sequences[header] = sequence
end

if not duplicates.empty? then
    duplicates.each { |k,v| puts "[WARN]: #{k} present #{v} times in input file"}
end

open(outfile, "w") do |outfh|
    outfh.write("#NEXUS\n")
    outfh.write("BEGIN TAXA;\n")
    outfh.write("\tDIMENSIONS NTAX = #{sequences.size};\n")
    outfh.write("\tTAXLABELS\n")
    if irregular_in_head == true then
        sequences.each_key { |k| outfh.write("\t'#{k}'\n")}
    else
        sequences.each_key { |k| outfh.write("\t#{k}\n")}
    end
    outfh.write(";\n")
    outfh.write("END;\n\n")
    outfh.write("BEGIN CHARACTERS;\n")
    outfh.write("\tDIMENSIONS NCHAR = #{length};\n")
    outfh.write("\tFORMAT\n")
    outfh.write("\t\tDATATYPE = DNA\n")
    outfh.write("\t\tGAP=-\n")
    outfh.write("\t\tMISSING=?\n")
    outfh.write("\t;\n\n")
    outfh.write("MATRIX\n")
    if irregular_in_head == true then
        sequences.each { |k,v| outfh.write("\t'#{k}'\t#{v}\n")}
    else
        sequences.each { |k,v| outfh.write("\t#{k}\t#{v}\n")}
    end
    outfh.write(";\nEND;\n")
end
