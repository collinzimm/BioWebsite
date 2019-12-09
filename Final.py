import pandas as pd
import os.path
from flask import Flask, render_template, request
import socket

# create pd df's
df82 = pd.read_feather("./82file.ftr", columns=None, use_threads=True)
df98 = pd.read_feather("./98file.ftr", columns=None, use_threads=True)
# clean for quicker search
genes82 = df82.query('feature == "gene"')
genes98 = df98.query('feature == "gene"')
trans82 = df82.query('feature == "transcript"')
trans98 = df98.query('feature == "transcript"')
exon82 = df82.query('feature == "exon"')
exon98 = df98.query('feature == "exon"')


def gene_annotation(genename):
    """
    Takes a gene name as only argument.
    returns a list containing year,
    chromosome, strand, start & end positions,
    and gene category for each release.
    var names: gList82, gList98
    """
    result82 = genes82.loc[genes82['gene_id'] == genename]
    g82 = result82.to_string(index=False).split('\n')
    geneList82 = g82[1].split()
    gList82 = ["82", geneList82[0], geneList82[6], geneList82[3], geneList82[4], geneList82[12]]

    result98 = genes98.loc[genes98['gene_id'] == genename]
    g98 = result98.to_string(index=False).split('\n')
    geneList98 = g98[1].split()
    gList98 = ["98",geneList98[0], geneList98[6], geneList98[3], geneList98[4], geneList98[12]]
    return gList82, gList98


def transcript_annotation(genename):
    """
    takes a gene as an argument.
    returns a list containing,
    year, start, end, exon number and intron number.
    var names: tList82, tList98
    """

    transcript = trans82.loc[trans82['gene_id'] == genename]
    x82 = transcript.to_string(header=False, index=False, index_names=False).split('\n')
    transList82 = x82[0].split()
    exon = exon82.loc[exon82['gene_id'] == genename]
    tList82 = ["82", transList82[3], transList82[4], exon.shape[0], exon.shape[0]-1]

    transcript = trans98.loc[trans98['gene_id'] == genename]
    x98 = transcript.to_string(header=False, index=False, index_names=False).split('\n')
    transList98 = x98[0].split()
    exon = exon98.loc[exon98['gene_id'] == genename]
    tList98 = ["98", transList98[3], transList98[4], exon.shape[0], exon.shape[0]-1]
    return tList82, tList98

def summary():
    """
    no args, returns total transcripts,
    total genes, total transcripts,
    gene and transcript types and numbers
    for each type, for each release
    var names: transNum, geneNum,
    biotypeGeneDict82, biotypeTransDict82,
    biotypeGeneDict98, biotypeTransDict98
    TODO: can probably just save in a file these outputs since they will always be the same
    """
    # print("\nSummary")
    transNum = trans82.shape[0] + trans98.shape[0]
    geneNum = df82.loc[df82['feature'] == 'gene'].shape[0] + df98.loc[df98['feature'] == 'gene'].shape[0]
    geneCats = df82['gene_biotype'].drop_duplicates()
    geneCats.append(df98['gene_biotype'].drop_duplicates())
    arr = geneCats.values  # this contains all categories
    # gene categories and relevant gene and transcript numbers
    # print("total transcripts =", transNum, "total genes =", geneNum)
    biotypeGeneDict82 = {}
    biotypeTransDict82 = {}
    biotypeGeneDict98 = {}
    biotypeTransDict98 = {}
    # add gene to dict with number
    for gene in arr:
        num1 = trans82.loc[trans82['gene_biotype'] == gene].shape[0]
        biotypeGeneDict82[gene] = num1
        num2 = genes82.loc[genes82['gene_biotype'] == gene].shape[0]
        biotypeTransDict82[gene] = num2
        num3 = trans98.loc[trans98['gene_biotype'] == gene].shape[0]
        biotypeGeneDict98[gene] = num3
        num4 = genes98.loc[genes98['gene_biotype'] == gene].shape[0]
        biotypeTransDict98[gene] = num4
    return geneNum, transNum, biotypeGeneDict82, biotypeGeneDict98, biotypeTransDict82, biotypeTransDict98


app = Flask(__name__)
@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':  # this block is only entered when the form is submitted
        geneName = request.form.get('gene name')
        geneNum, transNum, gene82, gene98, trans82, trans98 = summary()
        return render_template('results.html',
                               gene = gene_annotation(geneName),
                               trans = transcript_annotation(geneName),
                               geneNum = geneNum,
                               transNum = transNum,
                               gene82 = gene82,
                               gene98 = gene98,
                               trans82 = trans82,
                               trans98 = trans98,
                               geneName = geneName)
    return render_template('home.html')


if __name__ == '__main__':
    hostIp = socket.gethostbyname_ex(socket.gethostname())[-1]
    app.run(host=hostIp[1], port=5010)

# 65(not very efficient) lines of code
# ENSMUSG00000005610
