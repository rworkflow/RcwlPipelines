#!/usr/bin/env python

from optparse import OptionParser  
import gzip, yaml, re, sys, logging
import os.path

'''
'   Amie Radenbaugh - 02/09/2015
'   Jeltje van Baren - 06/31/2015
'   Program name: "vcfToMageTab.py"
'''

def get_read_fileHandler(aFilename):
    '''
    ' Open aFilename for reading and return
    ' the file handler.  The file can be
    ' gzipped or not.
    '''
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'rb')
    else:
        return open(aFilename,'r')


def get_write_fileHandler(aFilename):
    '''
    ' Open aFilename for writing and return
    ' the file handler.  The file can be
    ' gzipped or not.
    '''
    if aFilename.endswith('.gz'):
        return gzip.open(aFilename,'wb')
    else:
        return open(aFilename,'w')


def getConfig(inputYaml):
    """Extracts parameters from input yaml file, returns a single dict"""

    with open(inputYaml) as handle:
        configY = yaml.load(handle.read())
    handle.close()

    # FIXME This can be a command line input, but cannot always be parsed from the VCF header.
    # do we have all inputs?
    expected = ['expDesign', 'expDesignOntology', 'expDesignFactorName', 'expDesignFactorType', 'investigationTitle',
        'personLastName', 'personFirstName', 'personMidInitial', 'personEmail', 'personAddress', 
         'personAffiliation', 'personRole', 'pubMedId', 'pubAuthors', 'pubTitle', 'pubStatus', 
         'expDescription', 'protocolNames', 'protocolTypes', 'protocolDescriptions', 
         'protocolOntologies', 'protocolParameters', 'ontologyName', 'ontologyFile', 
         'ontologyVersion']
    found = set(configY.keys())
    if set(expected).difference(found):
        for mis in set(expected).difference(found):
            sys.stderr.write("ERROR, missing %s in yaml inputs\n" % mis)
        sys.exit(1)

    if set(found).difference(expected):
        for toomuch in set(found).difference(expected):
            sys.stderr.write("WARNING, %s in yaml inputs does not match to an input variable, ignoring...\n" % toomuch)
    return configY

def getSampleParams(inline):
##SAMPLE=<ID=NORMAL,SampleTCGABarcode=TCGA-BI-A0VR-10A-01D-A10S-08,SampleName=TCGA-BI-A0VR-10A-01D-A10S-08,Individual=TCGA-BI-A0VR,Description="Normal sample",Platform=Illumina,Source=dbGaP,Accession=.,softwareName=<muTect,CallIndelsPipeline>,softwareVer=<119,65>,softwareParam=<.>,File=TCGA-BI-A0VR-10A-01D-A10S-08,SampleUUID=e717369f-73df-46ca-a7f2-749b84b96616,MetadataResource=https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/xml/uuid/e717369f-73df-46ca-a7f2-749b84b96616>
    # first replace commas between nested <>, then split on comma
    line = inline
    matchObj = re.findall("<.*?>", line)
    for i in  matchObj:
        if re.search(',', i):
            commaToSC = re.sub(',', ';', i)
            inline = re.sub(i, commaToSC, inline)
    return inline.split(",")

def noneClean(v):
    if v is None:
        return ""
    return v

def createSDRF(vcfFile, sdrfFilename, archiveName, varDict, anIsDebug):
    """Create SDRF format file"""
    sdrfFileHandler = get_write_fileHandler(sdrfFilename)

    # output the header lines
    headerList = ["Extract Name", "Comment [TCGA Barcode]", "Comment [is tumor]", "Material Type", "Annotation REF", "Comment [TCGA Genome Reference]"]
    # some default ones that will be empty for us
    headerList += ["Protocol REF", "Parameter Value [Vendor]", "Parameter Value [Catalog Name]", "Parameter Value [Catalog Number]"]
    # protocol for bam file name
    headerList += ["Protocol REF", "Comment [Derived Data File REF]", "Comment [TCGA CGHub ID]", "Comment [TCGA Include for Analysis]"]

    # protocol somatic_variant_detection
    headerList += ["Protocol REF", "Derived Data File", "Comment [TCGA Spec Version]", "Comment [TCGA Include for Analysis]"]
    headerList += ["Comment [TCGA Data Type]", "Comment [TCGA Data Level]", "Comment [TCGA Archive Name]"]

    # add the header lines
    sdrfFileHandler.write("\t".join(headerList) + "\n")

    # hard-coded values
    emptyValue = "->"

    rnaFasta = "GRCh37-lite"	# FIXME: This doesn't seem appropriate
    dnaFasta = None
    dnaNormalBarcode = None
    dnaNormalUUID = None
    dnaNormalBamFilename = None
    dnaNormalCgHubId = None
    dnaTumorBarcode = None
    dnaTumorUUID = None
    dnaTumorBamFilename = None
    dnaTumorCgHubId = None
    rnaNormalBarcode = None
    rnaNormalUUID = None
    rnaNormalBamFilename = None
    rnaNormalCgHubId = None
    rnaTumorBarcode = None
    rnaTumorUUID = None
    rnaTumorBamFilename = None
    rnaTumorCgHubId = None

    # open the file
    vcfFileHandler = get_read_fileHandler(vcfFile)
    for line in vcfFileHandler:
        # strip the carriage return and newline characters
        line = line.rstrip("\r\n")

        if (anIsDebug):
            logging.debug("vcfLine: %s", line)

        # if it is an empty line, then just continue
        if (line.isspace()):
            continue;

        # we need to extract the tcga spec that was used
        elif (line.startswith("##tcgaversion")):
            ##tcgaversion=1.0

            (key, value) = line.split("=")
            tcgaSpecVersion = value
        # we need to extract info from the reference tag in the header
        elif (line.startswith("##reference")):
            ##reference=<ID=GRCh37,Source=file:/inside/depot/fa/Homo_sapiens_assembly19.fasta>
            ##assembly=file:/inside/depot/fa/Homo_sapiens_assembly19.fasta

            if (anIsDebug):
                logging.debug(line)

            if ("GRCh37-lite" in line):
                dnaFasta = "GRCh37-lite"
            elif ("GRCh37" in line):
                dnaFasta = "GRCh37"
            elif ("Homo_sapiens_assembly19" in line):
                dnaFasta = "Homo_sapiens_assembly19.fasta"
            elif ("hg19" in line):
                dnaFasta = "hg19"
            elif ("hg18" in line):
                dnaFasta = "hg18"
            elif ("NCBI36" in line):
                dnaFasta = "36.1"
            elif ("NCBI37" in line):
                dnaFasta = "37"
            else:
                dnaFasta = "NA"

        # we need to extract info from the SAMPLE tag in the header
        elif (line.startswith("##SAMPLE")):
            ##SAMPLE=<ID=DNA_NORMAL,SampleUUID=5218f2a6-5b5d-4e1a-aaf1-a0cd5e536571,SampleTCGABarcode=TCGA-IB-7646-10A-01D-2154-08,Individual=TCGA-IB-7646,
            #Description="Blood Derived Normal DNA",File="/inside/depot4/users/aradenba/data/hg19/paad/bams/wxs/TCGA-IB-7646-10A-01D-2154-08.bam",
            #Platform="Illumina HiSeq 2000",Source=CGHub,Accession=a6fd3469-b442-416d-aae9-991092fcb579,SequenceSource=WXS>
            ##SAMPLE=<ID=DNA_TUMOR,SampleUUID=e26b2473-9399-4eb8-9574-5b934b560740,SampleTCGABarcode=TCGA-IB-7646-01A-11D-2154-08,Individual=TCGA-IB-7646,
            #Description="Primary solid Tumor DNA",File="/inside/depot4/users/aradenba/data/hg19/paad/bams/wxs/TCGA-IB-7646-01A-11D-2154-08.bam",
            #Platform="Illumina HiSeq 2000",Source=CGHub,Accession=58287578-15b7-41d2-a5ac-70aaeac1fa40,SequenceSource=WXS>
            ##SAMPLE=<ID=RNA_TUMOR,SampleUUID=c3d463e2-8625-403f-af7b-63029bcb6eab,SampleTCGABarcode=TCGA-IB-7646-01A-11R-2156-07,Individual=TCGA-IB-7646,
            #Description="Primary solid Tumor RNA",File="/inside/depot4/users/aradenba/data/hg19/paad/bams/rna/TCGA-IB-7646-01A-11R-2156-07.bam",
            #Platform="Illumina HiSeq 2000",Source=CGHub,Accession=f4b705d6-29e2-4b83-9a99-3ef36d23d39b,SequenceSource=RNA-Seq>

            sampleLine = line[0:(len(line)-1)]
            sampleLine = sampleLine[len("##SAMPLE=<"):len(sampleLine)]
            sampleParamsList = getSampleParams(sampleLine)
            sampleType = -1
            sampleBarcode = ""
            sampleUUID = ""
            sampleBamFilename = ""
            sampleCgHubId = ""
            isRNA = False

            # create a dictionary of existing params
            for param in sampleParamsList:
                if ("software" in param):
                    continue;
                (key, value) = param.split("=")

                # remove the quotes
                if (value.startswith("\"")):
                    value = value[1:(len(value)-1)]

                # keep track of the barcode, uuid, and sequence source
                if (key == "SampleTCGABarcode"):
                    sampleBarcode = value
                    sampleBarcodeList = sampleBarcode.split("-")
                    sampleType = int(sampleBarcodeList[3][:2])
                    analyteType = sampleBarcodeList[4][2:]
                    if (analyteType == "R"):
                        isRNA = True
                elif (key == "SampleUUID"):
                    sampleUUID = value
                elif (key == "File"):
                    sampleBamFilename = value
                elif (key == "Accession"):
                    sampleCgHubId = value

            if (anIsDebug):
                logging.debug("SampleTCGABarcode=%s, SampleUUID=%s, SampleFile=%s, CGHubId=%s", sampleBarcode, sampleUUID, sampleBamFilename, sampleCgHubId)
            # the sample types from 0-9 are tumor and 10-19 are normal, 20 and above are control samples
            if (sampleType != -1):
                if (sampleType < 10):
                    if (isRNA):
                        rnaTumorBarcode = sampleBarcode
                        rnaTumorUUID = sampleUUID
                        rnaTumorBamFilename = sampleBamFilename
                        rnaTumorCgHubId = sampleCgHubId
                    else:
                        dnaTumorBarcode = sampleBarcode
                        dnaTumorUUID = sampleUUID
                        dnaTumorBamFilename = sampleBamFilename
                        dnaTumorCgHubId = sampleCgHubId
                elif (sampleType >= 10 and sampleType < 20):
                    if (isRNA):
                        rnaNormalBarcode = sampleBarcode
                        rnaNormalUUID = sampleUUID
                        rnaNormalBamFilename = sampleBamFilename
                        rnaNormalCgHubId = sampleCgHubId
                    else:
                        dnaNormalBarcode = sampleBarcode
                        dnaNormalUUID = sampleUUID
                        dnaNormalBamFilename = sampleBamFilename
                        dnaNormalCgHubId = sampleCgHubId
                else:
                    logging.critical("Traceback:  Unexpected sample type %s", sampleType)
                    sys.exit(1)

                if (anIsDebug):
                    logging.debug("normalBarcode=%s, normalUUID=%s, tumorBarcode=%s, tumorUUID=%s, rnaNormalBarcode=%s, rnaNormalUUID=%s, rnaTumorBarcode=%s, rnaTumorUUID=%s", dnaNormalBarcode, dnaNormalUUID, dnaTumorBarcode, dnaTumorUUID, rnaNormalBarcode, rnaNormalUUID, rnaTumorBarcode, rnaTumorUUID)

        # these are other header lines that we don't need, so just continue
        elif (line.startswith("#")):
            continue;

    # now we are to the data
    vcfFileHandler.close()
    protocolSomaticVariants = None
    protocolList = noneClean(varDict["protocolNames"]).split(",")
    for protocol in protocolList:
            protocolSomaticVariants = protocol

    vcfFilename = os.path.basename(vcfFile)
    # output one line per patient sample
    # normal and tumor DNA are required
    # check if RNA was available

    #################
    # Normal DNA
    #################
    # uuid, barcode, isTumor?, DNA/RNA, annotationRef, tcgaRef
    dnaNormalList = [dnaNormalUUID, dnaNormalBarcode, "no", "DNA", emptyValue, dnaFasta]
    # protocolRef:default, vendor, catalogName, CatalogNumber
    dnaNormalList += [emptyValue, emptyValue, emptyValue, emptyValue]
    # protocolRef:bamfile, bamfile, cgHubId, used?
    dnaNormalList += [emptyValue, dnaNormalBamFilename, dnaNormalCgHubId, "yes"]

    # protocol:somatic_variant, vcfFilename, tcgaSpecVersion, used?, dataType, dataLevel, archiveName
    dnaNormalList += [protocolSomaticVariants, vcfFilename, tcgaSpecVersion, "yes", "Mutations", "Level 2", archiveName]
    # write to file
    print("DNA Normal List", dnaNormalList)
    sdrfFileHandler.write("\t".join(dnaNormalList) + "\n")

    #################
    # Tumor DNA
    #################
    # uuid, barcode, isTumor?, DNA/RNA, annotationRef, tcgaRef
    dnaTumorList = [dnaTumorUUID, dnaTumorBarcode, "yes", "DNA", emptyValue, dnaFasta]
    # protocolRef:default, vendor, catalogName, CatalogNumber
    dnaTumorList += [emptyValue, emptyValue, emptyValue, emptyValue]
    # protocolRef:bamfile, bamfile, cgHubId, used?
    dnaTumorList += [emptyValue, dnaTumorBamFilename, dnaTumorCgHubId, "yes"]

    # protocol:somatic_variant, vcfFilename, tcgaSpecVersion, used?, dataType, dataLevel, archiveName
    dnaTumorList += [protocolSomaticVariants, vcfFilename, tcgaSpecVersion, "yes", "Mutations", "Level 2", archiveName]

    sdrfFileHandler.write("\t".join(dnaTumorList) + "\n")

    #################
    # Normal RNA
    #################
    if (rnaNormalBarcode != None):
        # uuid, barcode, isTumor?, DNA/RNA, annotationRef, tcgaRef
        rnaNormalList = [rnaNormalUUID, rnaNormalBarcode, "no", "RNA", emptyValue, rnaFasta]
        # protocolRef:default, vendor, catalogName, CatalogNumber
        rnaNormalList += [emptyValue, emptyValue, emptyValue, emptyValue]
        # protocolRef:bamfile, bamfile, cgHubId, used?
        rnaNormalList += [emptyValue, rnaNormalBamFilename, rnaNormalCgHubId, "yes"]

        # protocol:somatic_variant, vcfFilename, tcgaSpecVersion, used?, dataType, dataLevel, archiveName
        rnaNormalList += [protocolSomaticVariants, vcfFilename, tcgaSpecVersion, "yes", "Mutations", "Level 2", archiveName]

        sdrfFileHandler.write("\t".join(rnaNormalList) + "\n")

    #################
    # Tumor RNA
    #################
    if (rnaTumorBarcode != None):
        # uuid, barcode, isTumor?, DNA/RNA, annotationRef, tcgaRef
        rnaTumorList = [rnaTumorUUID, rnaTumorBarcode, "yes", "RNA", emptyValue, rnaFasta]
        # protocolRef:default, vendor, catalogName, CatalogNumber
        rnaTumorList += [emptyValue, emptyValue, emptyValue, emptyValue]
        # protocolRef:bamfile, bamfile, cgHubId, used?
        rnaTumorList += [emptyValue, rnaTumorBamFilename, rnaTumorCgHubId, "yes"]

        # protocol:somatic_variant, vcfFilename, tcgaSpecVersion, used?, dataType, dataLevel, archiveName
        rnaTumorList += [protocolSomaticVariants, vcfFilename, tcgaSpecVersion, "yes", "Mutations", "Level 2", archiveName]

        sdrfFileHandler.write("\t".join(rnaTumorList) + "\n")

    # close the file
    sdrfFileHandler.close()


def createIDF(idfFilename, sdrfFilename, varDict):
    """Create IDF format output. This looks very similar to the expected yaml config input, but contains the sdrf filename and a title"""
    idfFileHandler = get_write_fileHandler(idfFilename)

    # output the experimental design lines
    # FIXME: Title is hardcoded
    idfFileHandler.write("\t".join(["Investigation Title", noneClean(varDict["investigationTitle"])]) + "\n")
    idfFileHandler.write("\t".join(["Experimental Design", noneClean(varDict["expDesign"])]) + "\n")
    idfFileHandler.write("\t".join(["Experimental Design Term Source REF", noneClean(varDict["expDesignOntology"])]) + "\n")
    idfFileHandler.write("\t".join(["Experimental Factor Name", noneClean(varDict["expDesignFactorName"])]) + "\n")
    idfFileHandler.write("\t".join(["Experimental Factor Type", noneClean(varDict["expDesignFactorType"])]) + "\n")
    idfFileHandler.write("\n")

    # output the person lines
    idfFileHandler.write("\t".join(["Person Last Name", noneClean(varDict["personLastName"])]) + "\n")
    idfFileHandler.write("\t".join(["Person First Name", noneClean(varDict["personFirstName"])]) + "\n")
    idfFileHandler.write("\t".join(["Person Mid Initials", noneClean(varDict["personMidInitial"])]) + "\n")
    idfFileHandler.write("\t".join(["Person Email", noneClean(varDict["personEmail"])]) + "\n")
    idfFileHandler.write("\t".join(["Person Address", noneClean(varDict["personAddress"])]) + "\n")
    idfFileHandler.write("\t".join(["Person Affiliation", noneClean(varDict["personAffiliation"])]) + "\n")
    idfFileHandler.write("\t".join(["Person Roles", noneClean(varDict["personRole"])]) + "\n")
    idfFileHandler.write("\n")

    # output the publication lines
    idfFileHandler.write("\t".join(["PubMed ID", str(noneClean(varDict["pubMedId"]))]) + "\n")
    idfFileHandler.write("\t".join(["Publication Author List", noneClean(varDict["pubAuthors"])]) + "\n")
    idfFileHandler.write("\t".join(["Publication Title", noneClean(varDict["pubTitle"])]) + "\n")
    idfFileHandler.write("\t".join(["Publication Status", noneClean(varDict["pubStatus"])]) + "\n")
    idfFileHandler.write("\t".join(["Experiment Description", noneClean(varDict["expDescription"])]) + "\n")
    idfFileHandler.write("\n")

    # output the protocol lines
    idfFileHandler.write("\t".join(["Protocol Name", "\t".join(noneClean(varDict["protocolNames"]).split(","))]) + "\n")
    idfFileHandler.write("\t".join(["Protocol Type", "\t".join(noneClean(varDict["protocolTypes"]).split(","))]) + "\n")
    idfFileHandler.write("\t".join(["Protocol Description", "\t".join(noneClean(varDict["protocolDescriptions"]).split(","))]) + "\n")
    idfFileHandler.write("\t".join(["Protocol Term Source REF", "\t".join(noneClean(varDict["protocolOntologies"]).split(","))]) + "\n")
    idfFileHandler.write("\t".join(["Protocol Parameters", "\t".join(noneClean(varDict["protocolParameters"]).split(","))]) + "\n")
    idfFileHandler.write("\n")

    # output the sdrf line
    sdrfBasename = os.path.basename(sdrfFilename)
    idfFileHandler.write("\t".join(["SDRF Files", sdrfBasename]) + "\n")
    idfFileHandler.write("\n")

    # output the ontology lines
    idfFileHandler.write("\t".join(["Term Source Name", noneClean(varDict["ontologyName"])]) + "\n")
    idfFileHandler.write("\t".join(["Term Source File", noneClean(varDict["ontologyFile"])]) + "\n")
    idfFileHandler.write("\t".join(["Term Source Version", noneClean(varDict["ontologyVersion"])]) + "\n")

    # close the file
    idfFileHandler.close()



def main():

    # FIXME: This script does not yet work with the input xml
    # create the usage statement
    usage = """usage: python %prog <idf.yaml> <input.vcf> outputBase

Inputs are a yaml format file that contains standard parameters, used to create the IDF, 
an input vcf or just the header, and a base output filename, which will be used
to create outputBase.idf and outputBase.sdrf files
"""
    i_cmdLineParser = OptionParser(usage=usage)

    i_cmdLineParser.add_option("-l", "--log", dest="logLevel", default="WARNING", metavar="LOG", help="the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), %default by default")
    i_cmdLineParser.add_option("-g", "--logFilename", dest="logFilename", metavar="LOG_FILE", help="the name of the log file, STDOUT by default")

    if len(sys.argv) != 4:
        i_cmdLineParser.print_help()
        sys.exit(1)

    # get the required parameters
    (i_cmdLineOptions, i_cmdLineArgs) = i_cmdLineParser.parse_args()
    inputYaml = str(i_cmdLineArgs[0])
    inputVcf = str(i_cmdLineArgs[1])
    outBase = str(i_cmdLineArgs[2])
    idfFilename = outBase + '.idf'
    sdrfFilename = outBase + '.sdrf'
    archiveName = 'WhatShouldIBe'	# FIXME: this file is listed in the SDRF output, not sure what it should contain

    # parse the idf config
    args = getConfig(inputYaml)

    # get the optional params with default values
    i_logLevel = i_cmdLineOptions.logLevel
    i_logFilename = None
    # try to get any optional parameters with no defaults
    if (i_cmdLineOptions.logFilename != None):
        i_logFilename = str(i_cmdLineOptions.logFilename)

    # assuming loglevel is bound to the string value obtained from the
    # command line argument. Convert to upper case to allow the user to
    # specify --log=DEBUG or --log=debug
    i_numericLogLevel = getattr(logging, i_logLevel.upper(), None)
    if not isinstance(i_numericLogLevel, int):
        raise ValueError("Invalid log level: '%s' must be one of the following:  DEBUG, INFO, WARNING, ERROR, CRITICAL", i_logLevel)

    # set up the logging
    if (i_logFilename != None):
        logging.basicConfig(level=i_numericLogLevel, filename=i_logFilename, filemode='w', format='%(asctime)s\t%(levelname)s\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=i_numericLogLevel, format='%(asctime)s\t%(levelname)s\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # set the debug
    i_debug = (i_numericLogLevel < logging.WARNING)

    # do some debugging
    if (i_debug):
        logging.debug("idfFilename=%s", idfFilename)
        logging.debug("sdrfFilename=%s", sdrfFilename)
        logging.debug("i_logLevel=%s", i_logLevel)
        logging.debug("logFile=%s", i_logFilename)

    # create the output files
    createIDF(idfFilename, sdrfFilename, args)
    createSDRF(inputVcf, sdrfFilename, archiveName, args, i_debug)

    return

main()
sys.exit(0)
