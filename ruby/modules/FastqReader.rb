class FastqReader
    attr_writer :filename
    
    def initialize(filename)
        @filename = filename
        if not File.exists? filename then
            raise Errno::ENOENT
        end
    end
    
    def parse_records()
        header = '', sequence = '', quality = ''
        is_qual = false
        open(@filename).each do |line|
            line.chomp!
            next if line.empty?
            if line[/^@/] then
                if not header.empty? then
                    yield header, sequence, quality
                    sequence = '', quality = ''
                    is_qual = false
                end
                header = line[1..-1]
            elsif line =~ /^\+$/ then
                is_qual = true
            else
                if is_qual then
                    quality << line
                else
                    sequence << line
                end
            end
        end
        yield header, sequence, quality
    end
end