class FastaReader
    attr_accessor :filename
    
    def initialize(filename)
        @filename = filename
    end
    
    def parse_records()
        header = '', sequence = ''
        open(@filename).each do |line|
            line.chomp!
            next if line.empty?
            if line[/^>/] then
                if not header.empty? then
                    yield header, sequence if not sequence.empty?
                    sequence = ''
                end
                header = line[1..-1]
            else
                sequence << line
            end
        end
        yield header, sequence if not sequence.empty?
    end
end