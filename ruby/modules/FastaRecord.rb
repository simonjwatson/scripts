class FastaRecord

    attr_accessor :header, :sequence, :length

    def initialize(header, sequence)
        raise RuntimeError if sequence.empty?
        @header = header
        @sequence = sequence
        @length = sequence.length
    end
    
    def gc_perc
       freq = @sequence.count("gcGC") /Float(@length)
       return freq * 100
    end
    
    def reverse_complement
        rc = Hash["a" => "t", "c" => "g", "g" => "c", "t" => "a", "n" => "n"]
        rc_seq = ''
        @sequence.each_char { |base|
            rc_seq << rc[base.downcase]
        }
        return rc_seq.reverse
    end
end