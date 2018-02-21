from pylatex import Document, Section, Subsection, Command
from pylatex.utils import italic, NoEscape
import json
import os
import glob
import re


order_number = input("Which order to generate data for? ")
order_number = int(order_number)

author_set = set()
def fill_document(doc):
    """Add a section, a subsection and some text to the document.

    :param doc: the document
    :type doc: :class:`pylatex.document.Document` instance
    """
    for file in glob.glob("../../data/*/*.json"):
        with open(file,"r") as json_file:
            data = json.load(json_file)

        if data["info"]["order_number"] == order_number:
            with doc.create(Section(data["gene_name"])):
                doc.append(data["description"])
                with doc.create(Subsection('Author')):
                    doc.append(data["author"]["name"])
                    doc.append(',')
                    doc.append(data["author"]["email"])

if __name__ == '__main__':
    # Basic document
    document_title = ("Free Genes Twist Order %s" % order_number)
    doc = Document('basic')
    
    doc.preamble.append(Command('title', document_title))
    
    for file in glob.glob("../../data/*/*.json"):
        with open(file,"r") as json_file, open('fasta-documentation.txt',"a") as fasta_file:
            data = json.load(json_file)
        if data["info"]["order_number"] == order_number:
            author_set.add(data["author"]["name"])
            ### Write the fasta documentation
            fh = open('fasta-documentation.txt',"a")
            fh.write(">" + data["gene_name"] + "\n")
            fh.write(data["sequence"]["optimized_sequence"] + "\n") 
            fh.close()


    authors = ""
    for n in author_set:
        authors = authors + n + ", "
    print (authors)
    print (author_set)

    doc.preamble.append(Command('author', authors))
    doc.preamble.append(Command('date', NoEscape(r'\today')))
    doc.append(NoEscape(r'\maketitle'))
    fill_document(doc)
    doc.generate_tex()
    doc.generate_pdf('basic_maketitle', clean_tex=False)
    tex = doc.dumps()  # The document as string in LaTeX syntax
