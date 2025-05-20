# from dragibus.utils import get_sequence
from multiprocessing.dummy import Value
import pybedtools
import logging

class DragibusException(Exception):
    def __init__(self,key="Other dragibus exception",message=""):
        super().__init__(message)
        self.key = key

class WrongChromosomeExonError(DragibusException):
    def __init__(self,message="Trying to add new exon on different chromosome - skipping"):
        super().__init__("Exon with incorrect chromosome",message)

class WrongChromosomeTranscriptError(DragibusException):
    def __init__(self,message="Trying to add new transcript on different chromosome - skipping"):
        super().__init__("Transcript with incorrect chromosome",message)

class WrongStrandTranscriptError(DragibusException):
    def __init__(self,message="Trying to add new transcript with different strand - skipping"):
        super().__init__("Transcript with incorrect strand",message)

class WrongStrandExonError(DragibusException):
    def __init__(self,message="Trying to add new exon with different strand - skipping"):
        super().__init__("Exon with incorrect strand",message)

class WrongExonNumbering(DragibusException):
    def __init__(self,message="exon_number attribute does not match with the one inferred by Dragibus - renumbering"):
        super().__init__("Incorrect exon numbering",message)
    pass

class InvalidCoordinatesError(DragibusException):
    # Coordinates are not int
    def __init__(self):
        super().__init__("Entry with invalid coordinates","Invalid coordinates - skipping feature")

class InvalidStrandError(DragibusException):

    pass

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
    return "; ".join([f'{k} "{str(v)}"' for k, v in attributes.items()]) + ";"

class Feature:

    def __init__(self, fixed_fields, attributes):
        self.chr = fixed_fields[0]
        self.source = fixed_fields[1]
        self.type = fixed_fields[2]
        try:
            self.start = int(fixed_fields[3])
            self.end = int(fixed_fields[4])
        except ValueError:
            raise(InvalidCoordinatesError)
        self.score = fixed_fields[5]
        self.strand = fixed_fields[6]
        self.frame = fixed_fields[7]
        self.attributes = attributes
        self.length = self.end - self.start + 1

        try:
            if self.strand != '+' and self.strand != '-':
                raise(InvalidStrandError)
        except InvalidCoordinatesError as e:
            raise e
        except InvalidStrandError as e:
            raise e
    
    def format(self):
        s = ("\t".join([
            self.chr,
            self.source,
            self.type,
            str(self.start),
            str(self.end),
            str(self.score),
            self.strand,
            self.frame,
            attributes_to_str(self.attributes)
        ]))
        s+="\n"
        return(s)

class Gene(Feature):

    def __init__(self, fixed_fields, attributes):
        if isinstance(attributes,str):
            attributes = read_attributes(attributes)
        else:
            assert isinstance(attributes,dict)

        try:
            super(Gene, self).__init__(fixed_fields,attributes)
        except Exception as e:
            raise(e)
        self.id = self.attributes['gene_id']
        self.transcripts = []
    
    def add_transcript(self,t):
        self.transcripts.append(t)

        if t.chr != self.chr:
            raise  WrongChromosomeTranscriptError()
        if t.strand != self.strand:
            raise WrongStrandTranscriptError()

        self.start = min(t.start,self.start)
        self.end = max(t.end,self.end)

class Transcript(Feature):

    def __init__(self, fixed_fields, attributes):
        if isinstance(attributes,str):
            attributes = read_attributes(attributes)
        else:
            assert isinstance(attributes,dict)

        try:
            super(Transcript, self).__init__(fixed_fields,attributes)
        except Exception as e:
            raise(e)

        self.id = self.attributes['transcript_id']
        self.gene_id = self.attributes['gene_id']
        self.exons = []
        self.introns = []
        self.cdna_length = None
        self.transcript_length = None
        self.polyA = None
        self.CDS = False
        self.attributes["n_exons"] = 0

    def compute_length(self):
        self.transcript_length = self.end - self.start + 1 # to check

    def add_exon(self,e):
        # assert e.start >= self.end or e.end <= self.start  # Check that the transcript does not already covers the exons
        if e.chr != self.chr:
            raise WrongChromosomeExonError()
        if e.strand != self.strand:
            raise WrongStrandExonError()

        self.exons.append(e)
        self.start = min(self.start,e.start)
        self.end = min(self.end,e.end)
        self.attributes["n_exons"] += 1

    def sort_exons(self):
        wrong_numbering = False
        if len(self.exons) > 0:
            self.exons.sort(key=lambda x:x.start)
            self.start = self.exons[0].start
            self.end = self.exons[-1].end

            # Add exon number attribute if it's not present in the gtf
            for i,e in enumerate(self.exons):
                if not hasattr(e,'exon_number'):
                    e.exon_number = i+1
                else:
                    # Raise exception is dragibus has a different exon numbering
                    # Can happen especially when exons where ignored of of incorrect strand or other reason
                    if e.exon_number != i+1:
                        wrong_numbering = True
                        e.exon_number = i+1
                                       
            # Mark internal exons
            for i,e in enumerate(self.exons):
                if i==0 or i==(len(self.exons)-1):
                    e.is_internal = False
                else:
                    e.is_internal = True

        if wrong_numbering:
            raise WrongExonNumbering

    def sort_introns(self):
        self.introns.sort(key=lambda x:x.start)

    def add_introns(self):
        for i in range(0,len(self.exons)-1):
            start = self.exons[i].end + 1
            end = self.exons[i+1].start - 1

            attributes = {k:v for (k,v) in self.attributes.items() if k != 'exon_number'}

            try:
                f = Intron([self.chr,self.source,"intron",start,end,self.score,self.strand,self.frame],attributes)
            except Exception as e:
                raise(e)

            self.introns.append(f)

    def compute_cdna_length(self):
        if not self.cdna_length:
            self.cdna_length = sum([e.end - e.start + 1 for e in self.exons])

class Exon(Feature):

    def __init__(self, fixed_fields, attributes):
        if isinstance(attributes,str):
            attributes = read_attributes(attributes)
        else:
            assert isinstance(attributes,dict)

        try:
            super(Exon, self).__init__(fixed_fields,attributes)
        except Exception as e:
            raise e
        
        if 'exon_number' in self.attributes.keys():
            self.exon_number =  int(self.attributes['exon_number'])
        self.gene_id = self.attributes['gene_id']
        self.transcript_id = self.attributes['transcript_id']

class Intron(Feature):
    def __init__(self, fixed_fields, attributes):
        if isinstance(attributes,str):
            attributes = read_attributes(attributes)
        else:
            assert isinstance(attributes,dict)

        try:
            super(Intron, self).__init__(fixed_fields,attributes)
        except Exception as e:
            raise e
        
        self.id = "_".join([self.chr,str(self.start),str(self.end),str(self.strand)])
        self.canonic = None

    def get_donor_acceptor_coords(self):
        # TODO Check chromosome boundaries
        if self.strand =="+":
            return(((self.start-1,self.start+1),(self.end-2,self.end)))
        else:
            return(((self.end-2,self.end),(self.start-1,self.start + 1 )))

class CDS(Feature):
    def __init__(self, fixed_fields, attributes):
        if isinstance(attributes,str):
            attributes = read_attributes(attributes)
        else:
            assert isinstance(attributes,dict)

        try:
            super(CDS, self).__init__(fixed_fields,attributes)
        except Exception as e:
            raise e
        
        self.transcript_id = self.attributes['transcript_id']

def find_canonic_introns(transcripts,fasta):
    nb_introns = 0
    dimer_coords = set()

    with(open("canonic_introns_dragibus.txt",'w')) as outf:
        with(open("noncanonic_introns_dragibus.txt",'w')) as outf2:
            # Collect the coordinates of donor / acceptor sequences for each intron
            for t in transcripts.values():
                # print(t.strand)
                # print(t.id)
                for i in t.introns:
                    # print("i :"+i.strand)                    
                    
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

                    # if i.id=="12_5429406_5429503":
                    #     print((i.chr,coords[0][0],coords[0][1],i_id + "_donor",0,i.strand))
                    #     print((i.chr,coords[1][0],coords[1][1],i_id + "_acceptor",0,i.strand))

            # Get sequence for each dimer using bedtools
            bed = pybedtools.BedTool(list(dimer_coords))
            # print(bed)
            fa = bed.sequence(fi=fasta,s=True,name=True).seqfn

            donor_seq = dict() # dict intron_id -> str
            acceptor_seq = dict()

            with open(fa) as res:
                for line in res:
                    if line.startswith(">"):
                        header = line[1:].split("::")[0]
                        id = "_".join(header.split("_")[:-1])
                        if header.split("_")[-1].startswith('donor'):
                            dir = "donor"
                        elif header.split("_")[-1].startswith('acceptor'):
                            dir = "acceptor"
                    else:
                        if dir == "donor":
                            donor_seq[id] = line.strip().upper()
                        elif dir == "acceptor":
                            acceptor_seq[id] = line.strip().upper()
                        id = dir = ""
            # print(len(donor_seq))
            # Establish canonicity
            a=0
            for t in transcripts.values():
                for i in t.introns:
                    a += 1

                    if (acceptor_seq[i.id],donor_seq[i.id]) in (("AG","GT"),("AG","GC"),("AC","AT")):
                        i.canonic=True
                        i.attributes["is_canonical"] = "True"
                        outf.write(" ".join([i.chr,str(i.start),str(i.end),'"'+t.id+'";'])+"\n")
                    else:
                        i.canonic=False
                        i.attributes["is_canonical"] = "False"
                        outf2.write(" ".join([i.chr,str(i.start),str(i.end),'"'+t.id+'";'])+"\n")
                
                # Annotate if transcript has all introns canonical
                if {i.attributes['is_canonical'] for i in t.introns} == {"True"}:
                    t.attributes["all_introns_canonical"] = "True"
                else:
                    t.attributes["all_introns_canonical"] = "False"




