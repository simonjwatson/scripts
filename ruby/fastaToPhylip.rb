#! /usr/bin/ruby

require_relative 'modules/FastaReader'

def edit_header(header)
	return header.gsub(/\s/, '_')
end


if ARGV.length != 2 then
	puts "[USAGE]: fastaToPhylip.rb <infile.fas> <outfile.phy>"
	Process.exit
end

infile = ARGV[0]
outfile = ARGV[1]
samples = {}
len = 0

fasta = FastaReader.new(infile)
fasta.parse_records do |header, sequence|
	samples[header] = sequence
	l = sequence.length
	if len == 0 then
		len = l
	elsif l != len
		puts "[ERROR]: Mismatch in sequence lengths in FASTA"
		Process.exit
	end
end

open(outfile, "w") do |outfh|
	outfh.puts("#{samples.keys.length} #{len}")
	samples.each { |key, value| outfh.puts("#{edit_header(key)}\t#{value}") }
end
