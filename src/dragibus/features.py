# from dragibus.utils import get_sequence
import pybedtools

def read_attributes(attributes, quote_char='\"',missing_value=""):
    # potential tip to reduce memory usage using interning : 
    # https://github.com/openvax/gtfparse/blob/18cc0fde7498ad72c87610268f60826fd89e442b/gtfparse/attribute_parsing.py#L57
    
    if isinstance(attributes,str):
    
        attributes_dict = dict()

        for field in attributes.split(";"):
            field = field.strip().split(" ")

            if len(field) != 2:
                continue
            
            name,value = field
            value = value.replace(quote_char, "") if value.startswith(quote_char) else value
            attributes_dict[name] = value

        return(attributes_dict)
    elif isinstance(attributes,dict):
        return(attributes)
    else:
        raise Exception()

def attributes_to_str(attributes):
    return("; ".join([" ".join([k,'"'+v+'"']) for k,v in attributes.items()])+";")

class Feature:

    def __init__(self, fixed_fields, attributes):
        self.chr = fixed_fields[0]
        self.source = fixed_fields[1]
        self.type = fixed_fields[2]
        self.start = int(fixed_fields[3])
        self.end = int(fixed_fields[4])
        self.score = fixed_fields[5]
        self.strand = fixed_fields[6]
        self.frame = fixed_fields[7]
        self.attributes = attributes
        self.length = self.end - self.start + 1

class Gene(Feature):

    def __init__(self, fixed_fields, attributes):
        if isinstance(attributes,str):
            attributes = read_attributes(attributes)
        else:
            assert isinstance(attributes,dict)

        super(Gene, self).__init__(fixed_fields,attributes)
        self.id = self.attributes['gene_id']
        self.transcripts = []
    
    def add_transcript(self,t):
        self.transcripts.append(t)
        assert t.chr == self.chr
        assert t.strand == self.strand

        self.start = min(t.start,self.start)
        self.end = max(t.end,self.end)

class Transcript(Feature):

    def __init__(self, fixed_fields, attributes):
        if isinstance(attributes,str):
            attributes = read_attributes(attributes)
        else:
            assert isinstance(attributes,dict)

        super(Transcript, self).__init__(fixed_fields,attributes)
        self.id = self.attributes['transcript_id']
        self.gene_id = self.attributes['gene_id']
        self.exons = []
        self.introns = []
        self.cdna_length = None
        self.transcript_length = None
        self.polyA = None

    def compute_length(self):
        self.transcript_length = self.end - self.start + 1 # to check

    def add_exon(self,e):
        self.exons.append(e)
        assert e.chr == self.chr
        assert e.strand == self.strand
        # assert e.start >= self.end or e.end <= self.start  # Check that the transcript does not already covers the exons
        self.start = min(self.start,e.start)
        self.end = min(self.end,e.end)

    def sort_exons(self):
        if len(self.exons) > 0:
            self.exons.sort(key=lambda x:x.start)
            self.start = self.exons[0].start
            self.end = self.exons[-1].end

            # Mark internal exons
            for i,e in enumerate(self.exons):
                if i==0 or i==(len(self.exons)-1):
                    e.is_internal = False
                else:
                    e.is_internal = True

    def sort_introns(self):
        self.introns.sort(key=lambda x:x.start)

    def add_introns(self):
        for i in range(1,len(self.exons)-1):
            start = self.exons[i].end + 1
            end = self.exons[i+1].start - 1

            attributes = {k:v for (k,v) in self.attributes.items() if k != 'exon_number'}

            f = Intron([self.chr,self.source,"intron",start,end,self.score,self.strand,self.frame],attributes)

            self.introns.append(f)

    def compute_cdna_length(self):
        if not self.cdna_length:
            self.cdna_length = sum({e.end - e.start for e in self.exons})

class Exon(Feature):

    def __init__(self, fixed_fields, attributes):
        if isinstance(attributes,str):
            attributes = read_attributes(attributes)
        else:
            assert isinstance(attributes,dict)

        super(Exon, self).__init__(fixed_fields,attributes)
        self.exon_number =  self.attributes['exon_number']
        self.gene_id = self.attributes['gene_id']
        self.transcript_id = self.attributes['transcript_id']

class Intron(Feature):
    def __init__(self, fixed_fields, attributes):
        if isinstance(attributes,str):
            attributes = read_attributes(attributes)
        else:
            assert isinstance(attributes,dict)

        super(Intron, self).__init__(fixed_fields,attributes)
        
        self.id = "_".join([self.chr,str(self.start),str(self.end)])
        self.canonic = None

    def get_donor_acceptor_coords(self):
        # TODO Check chromosome boundaries
        if self.strand =="+":
            return(((self.start-1,self.start+1),(self.end-2,self.end)))
        else:
            return(((self.end-2,self.end),(self.start-1,self.start + 1 )))



def find_canonic_introns(transcripts,fasta):
    nb_introns = 0
    dimer_coords = set()

    # Collect the coordinates of donor / acceptor sequences for each intron
    for t in transcripts.values():
        for i in t.introns:
            
            nb_introns += 1
            i_id = i.id
            # print(nb_introns)
            coords = i.get_donor_acceptor_coords()
            # dimer_coords.add((i.chr,coords[0][0],coords[0][1],i_id + "_donor",0,i.strand))
            # dimer_coords.add((i.chr,coords[1][0],coords[1][1],i_id + "_acceptor",0,i.strand))
            # if i.id == "2_67875_68406":
            # print((i.chr,coords[0][0],coords[0][1],i_id + "_donor",0,i.strand))
            dimer_coords.add((i.chr,coords[0][0],coords[0][1],i_id + "_donor",0,i.strand))
            dimer_coords.add((i.chr,coords[1][0],coords[1][1],i_id + "_acceptor",0,i.strand))

    # Get sequence for each dimer using bedtools
    bed = pybedtools.BedTool(list(dimer_coords))
    fa = bed.sequence(fi=fasta,s=True,name=True).seqfn

    donor_seq = dict() # dict intron_id -> str
    acceptor_seq = dict()

    with open(fa) as res:
        for line in res:
            if line.startswith(">"):
                header = line[1:].split("::")[0]
                id = "_".join(header.split("_")[0:3])
                if header.split("_")[3].startswith('donor'):
                    dir = "donor"
                elif header.split("_")[3].startswith('acceptor'):
                    dir = "acceptor"
            else:
                if dir == "donor":
                    donor_seq[id] = line.strip().upper()
                elif dir == "acceptor":
                    acceptor_seq[id] = line.strip().upper()
                id = dir = ""
    # Establish canonicity
    a=0
    for t in transcripts.values():
        for i in t.introns:
            a += 1 
            if (acceptor_seq[i.id],donor_seq[i.id]) in (("AG","GT"),("AG","GC"),("AC","AT")):
                i.canonic=True
            else:
                i.canonic=False





