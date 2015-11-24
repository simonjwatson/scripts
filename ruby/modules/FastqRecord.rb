class FastqRecord
    attr_accessor :header, :sequence, :quality
    
    def initialize(header, sequence, quality)
        raise RuntimeError => err if sequence.length != quality.length
        @header = header
        @sequence = sequence
        @quality = quality
    end
end
