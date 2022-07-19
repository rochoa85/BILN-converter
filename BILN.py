"""Script to manipulate BILN and convert to HELM (and viceversa)

From publication: BILN – A Human-readable Line Notation for Complex Peptides
Journal of Chemical Information and Modelling, 2022
"""

########################################################################################
# Authorship
########################################################################################

__credits__ = ["Thomas Fox", "Michael Bieler", "Peter Haebel", "Rodrigo Ochoa", "Stefan Peters", "Alexander Weber"]
__license__ = "MIT"
__version__ = "1.0"

########################################################################################
# Modules
########################################################################################

import copy
import re
import logging
import argparse
import sys

##########################################################################
# Main class
##########################################################################
class BILN:

    def __init__(self, biln=None, helm=None):
        """Construct a BILN or HELM molecule. They can be read from a file

        :param biln: The BILN format that will be converted to HELM
        :type biln: str
        :param helm: The HELM format that will be converted to BILN
        :type helm: str
        """

        self.polymerinfo = {"chains": [], "bonds": []}

        if biln is not None and helm is not None:
            logger and logger.error("Either biln or helm or none are allowed as input. Exit the script")
            exit(1)
        elif biln is not None:
            self.evalBILN(biln=biln)
        elif helm is not None:
            self.evalHELM(helm=helm)
        
    ############################################################
    @staticmethod
    def __removeBrackets(sequence):
        """Remove brackets around individual residues

        Args:
            :param sequence: the peptide sequence
            :type sequence: list

        Returns:
            :param newseq: new list without the brackets
            :type newseq: list
        """
        newseq = []
        for r in sequence:
            p = re.sub(r'\[(.*)\]', r'\1', r)  # remove brackets if necessary
            newseq.append(p)

        return newseq

    ############################################################
    @staticmethod
    def __splitHELM(helm):
        """split HELM string into parts

        Args:
            :param helm: the helm peptide sequence
            :type sequence: string

        Returns:
            :param helm_parts: a list containing relevant parts from helm
            :type helm_parts: list
        """

        helm_parts = []
        while len(helm):
            m = re.search('\$[^,;]', helm)  # this excludes any $ signs in the cxsmiles
            if m is None:
                m = re.search('\$', helm)  # search for all other $ signs

            if m is None:
                helm_parts.append(helm)  # done, no more parts
                break

            else:
                helm_parts.append(helm[:m.span()[0]])

                helm = helm[m.span()[0] + 1:]  # remove the current part from helm

        if len(helm_parts) == 4:
            helm_parts.append('')

        # helm part 1 == list of simple polymers - split further
        listOfSimplePolymers = helm_parts[0]
        if '|' in listOfSimplePolymers:
            listOfSimplePolymers = listOfSimplePolymers.split("|")
        else:
            listOfSimplePolymers = [listOfSimplePolymers]
        helm_parts[0] = listOfSimplePolymers

        # split the bond information into individual chunks
        listOfConnections = helm_parts[1]
        if listOfConnections != '':
            if '|' in listOfConnections:
                listOfConnections = listOfConnections.split("|")
            else:
                listOfConnections = [listOfConnections]
        else:
            listOfConnections = []

        helm_parts[1] = listOfConnections

        return helm_parts

    ############################################################
    def __toBILN(self):
        """generate BILN from polymerInfo

        Returns:
            :param biln: the generated biln from polymer dictionary
            :type biln: string
        """

        chains = copy.deepcopy(self.polymerinfo["chains"])
        bonds = copy.deepcopy(self.polymerinfo["bonds"])

        # move bond info into monomers
        for ibond, bond in enumerate(bonds):
            c1, res1, Rgroup1, c2, res2, Rgroup2 = bond
            chains[c1][res1] = "{}({},{})".format(chains[c1][res1], ibond + 1, Rgroup1)
            chains[c2][res2] = "{}({},{})".format(chains[c2][res2], ibond + 1, Rgroup2)

        # chain everything together
        listOfsimplePolymers = []
        for c in chains:
            poly = "-".join(c)
            listOfsimplePolymers.append(poly)

        biln = ".".join(listOfsimplePolymers)

        if len(biln) == 0:
            biln = None

        return biln

    ############################################################
    def __toHELM(self):
        """generate HELM from polymerInfo

        Returns:
            :param helm: the generated helm v2.0 from polymer dictionary
            :type helm: string
        """

        chains = copy.deepcopy(self.polymerinfo["chains"])
        bonds = copy.deepcopy(self.polymerinfo["bonds"])

        # brackets around residues if necessary
        for ic, c in enumerate(chains):
            for ir, res in enumerate(c):
                if len(res) > 1:
                    chains[ic][ir] = "[{}]".format(res)

        # generate HELM elements
        # compile chains
        listOfSimplePolymers = []
        for ic, c in enumerate(chains):
            poly = ".".join(c)
            listOfSimplePolymers.append('PEPTIDE{}{{{}}}'.format(ic + 1, poly))

        # compile bondinfo
        listOfConnections = []
        for bond in bonds:
            c1, r1, g1, c2, r2, g2 = bond
            bondInfo = 'PEPTIDE{},PEPTIDE{},{}:R{}-{}:R{}'.format(c1 + 1, c2 + 1, r1 + 1, g1, r2 + 1, g2)
            listOfConnections.append(bondInfo)

        listOfConnections = '|'.join(listOfConnections)
        listOfSimplePolymers = '|'.join(listOfSimplePolymers)

        if not listOfSimplePolymers:
            helm = None
        else:
            helm = "{}${}$$$V2.0".format(listOfSimplePolymers, listOfConnections)

        return helm

    ############################################################
    def evalHELM(self, helm):
        """Convert HELM string to polymerinfo

        This method fills the polymerinfo dictionary:
        polymerinfo = {"chains": [chain1, chain2, ...],
                       "bonds": [b1, b2, b3]}
        where chains contains the monomer abbreviation for each chain
        and bonds contains the "extra" bonds explicitly given in HELM or BILN as a list

        Returns:
            :param helm: the peptide helm string
            :type helm: string
        """

        # split the helm string into its individual parts
        try:
            listOfSimplePolymers, listOfConnections, groups, annotations, version = self.__splitHELM(helm)

        except (ValueError, IndexError):
            logger and logger.error('problem with HELM string - not enough sections: '.format(helm))
            logger and logger.error('need 5, have {}'.format(len(self.__splitHELM(helm))))

            return None

        if len(listOfSimplePolymers) == 0:
            logger and logger.error('No simple polymers in HELM string {}'.format(helm))
            return None

        # go through the list of simple polymers (first component of HELM string), parse them and put them into
        # polymerinfo["chains"]

        pattern = re.compile(r'{.*}')

        id = []
        polymer = []

        for idx, chain in enumerate(listOfSimplePolymers):

            # remove any whitespace
            chain = chain.strip()
            # split each polymer into name/identifier and sequence
            m = pattern.search(chain)
            if m is None:
                logger and logger.error('no sequence information found in simple polymer - check HELM')
                logger and logger.error('input: {}'.format(chain))
                return None

            # split each polymer into name/identifier and sequence
            s = m.span()
            ic = chain[:s[0]]  # identifier
            if ic[0:4] == 'CHEM':
                logger and logger.error(
                    'polymer contains an explicit chemical entity - missing or not recognized monomer:')
                logger and logger.error(chain)
                return None
            elif ic[0:7] != 'PEPTIDE':
                logger and logger.error('non-peptide chains in HELM - probably missing monomer')
                logger and logger.error('found: {}'.format(ic))
                return None
            else:
                ic = int(re.sub('PEPTIDE', '', ic))

            p = chain[s[0] + 1:s[1] - 1]  # sequence

            if not len(p):
                logger and logger.error('simple polymer {} is declared but not defined - has no length'.format(p))
                return None

            p = splitOutside(p, '.', '[]')  # split into individual monomers.

            p = self.__removeBrackets(p)  # remove brackets around monomer abbreviations with more than 1 letter

            id.append(ic)
            polymer.append(p)

        self.polymerinfo["chains"] = polymer

        # parse the bond information

        # bond information is of the type: 'PEPTIDE1,PEPTIDE2,1:R1-4:R3'
        # split into individual parts 'id1, id2, res1:rgroup1-res2:rgroup2'
        bonds = []
        if len(listOfConnections) > 0:
            for idx, conn in enumerate(listOfConnections):
                id1, id2, bond = conn.split(',')

                res1, rgroup1, res2, rgroup2 = re.split('[-:]', bond)

                # need some reformatting
                id1 = int(id1.replace('PEPTIDE', ''))
                id2 = int(id2.replace('PEPTIDE', ''))

                # translate back to current order of chains in polymerinfo["chains"]
                id1 = id.index(id1)
                id2 = id.index(id2)

                # start counting of residues at 0
                res1 = int(res1) - 1
                res2 = int(res2) - 1

                # get the number of the Rgroup (keep numbering 1-3)
                rgroup1 = int(rgroup1.replace('R', ''))
                rgroup2 = int(rgroup2.replace('R', ''))

                logger and logger.debug('bond :{}'.format(', '.join([str(id1), str(res1), str(rgroup1), str(id2), str(res2), str(rgroup2)])))

                bonds.append([id1, res1, rgroup1, id2, res2, rgroup2])

        self.polymerinfo["bonds"] = bonds

    ############################################################
    def evalBILN(self, biln):
        """Generate polymerinfo dictionary from BILN string

        This method fills the polymerinfo dictionary:
        polymerinfo = {"chains": [chain1, chain2, ...],
                       "bonds": [b1, b2, b3]}
        where chains contains the monomer abbreviation for each chain
        and bonds contains the "extra" bonds explicitly given in HELM or BILN as a list

        Returns:
            :param biln: the peptide biln string
            :type biln: string
        """

        # split BILN into chains
        chains = biln.split(".")

        bondinfo = {}
        listOfSimplePolymers = []

        for cidx, c in enumerate(chains):
            residues = c.split("-")

            # go through residues and collect the bond information needed
            polymer = []
            for ridx, r in enumerate(residues):
                match = re.findall(r"\((\d+),(\d+)\)", r)
                resname = re.sub(r"\(.*\)", "", r)
                if match:
                    for m in match:
                        bidx = m[0]
                        gidx = int(m[1])
                        try:
                            bondinfo[bidx].append(cidx)
                        except:
                            bondinfo[bidx] = [cidx]

                        bondinfo[bidx].append(ridx)
                        bondinfo[bidx].append(gidx)

                polymer.append(resname)
            listOfSimplePolymers.append(polymer)

        self.polymerinfo["chains"] = listOfSimplePolymers

        # compile bondinfo
        listOfConnections = []
        for key in bondinfo:
            if len(bondinfo[key]) != 6:
                logger and logger.error('error in bond information of BILN: bond {} not correct'.format(key))
                return None

            c1, r1, g1, c2, r2, g2 = bondinfo[key]

            bond = [c1, r1, g1, c2, r2, g2]
            listOfConnections.append(bond)

        self.polymerinfo["bonds"] = listOfConnections

    ############################################################
    def getHELM(self):
        """return HELM string from polymerinfo"""

        helm = self.__toHELM()

        return helm

    ############################################################
    def getBILN(self):
        """return BILN string from polymerinfo"""
        biln = self.__toBILN()

        return biln

##########################################################################
# Additional function
##########################################################################
def splitOutside(string, by, outside, keepMarker=True):
    """splits a string by delimiter only if outside of a given delimiter (bracket, quote, ...)
    The function keeps track of opening and closing parens, brackets and braces, as well as single and double quotes,
    and performs a replacement only outside of such bracketed and quoted substrings.
    Then replaces the non-bracketed/quoted search characters with another character which definitely doesn't appear
    in the string (e.g. ASCII/Unicode group-separator: chr(29) code), then do a simple string.split on that character.

    Args:
        :param string: string to be split
        :type string: str
        :param by: delimiter(s) by which to be split
        :type by: str
        :param outside: only split if outside of this
        :type outside: str
        :param keepMarker: if True keep the chunk marker, remove otherwise
        :type keepMarker: bool

    Returns:
        :param splitChains: split string as list
        :type splitChains: list
    """
    # by can be more than 1 character
    by = list(by)

    # if outside is only one character (e.g. ', "), double it for start and end
    if len(outside) == 1:
        outside = outside + outside

    # Special character
    GRPSEP = chr(29)

    out = ''
    inside = False
    for i in string:
        if i == outside[0]:
            if inside:
                if keepMarker: j = i
                else: j = ''
                inside = False
            else:
                inside = True
                if keepMarker: j = i
                else: j = ''
        elif i == outside[1]:
            inside = False
            if keepMarker: j = i
            else: j = ''
        else:
            if not inside and i in by: j = GRPSEP
            else: j = i
        out = out + j

    # Do the final split
    splitChains = out.split(GRPSEP)
    return splitChains

##########################################################################
def getInputsParser():
    """Constructs parser for inputs."""

    parser = argparse.ArgumentParser(add_help=False)

    inputTypeGroup = parser.add_mutually_exclusive_group(required=True)
    inputTypeGroup.add_argument(
        '--biln', type=str, metavar='text',
        required=False,
        help="BILN string to convert into HELM format.")
    inputTypeGroup.add_argument(
        '--helm', type=str, metavar='text',
        required=False,
        help="HELM string to convert into BILN format.")
    inputTypeGroup.add_argument(
        '--table_biln', type=str, metavar='filename',
        required=False,
        help="Table containing BILN strings to convert into HELM.")
    inputTypeGroup.add_argument(
        '--table_helm', type=str, metavar='filename',
        required=False,
        help="Table containing HELM strings to convert into BILN.")

    # A repeated-use log option parser.
    logOptions = parser.add_argument_group('Logging options')
    logOptions.add_argument(
        '--logfile', type=str, metavar='filename',
        required=False, default=None,
        help="Output messages to given logfile, default is stderr.")
    logOptions.add_argument(
        "-v", "--verbose", action="store_true",
        help="Increase output verbosity")

    return parser

##########################################################################
# Main function
##########################################################################
if __name__ == "__main__":
    # Read arguments
    useParser = argparse.ArgumentParser(
        description="""Conver BILN to HELM and viceversa""",
        parents=(getInputsParser(),
        ))
    args = useParser.parse_args()

    # Setup typical logger for messages to stderr.
    logStream = logging.StreamHandler(
        open(args.logfile, "w") if args.logfile else sys.stderr)
    logStream.setFormatter(logging.Formatter(
        "%(asctime)-10s %(levelname)s:%(message)s",
        datefmt="%H:%M:%S"))
    logger = logging.Logger(name="Program")
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    logger.addHandler(logStream)
    logger.debug("Invocation arguments: {}".format(args))

    # Options pipelines
    # 1. A single BILN input directly in the command line
    if args.biln:
        b = BILN(biln=args.biln)
        helm = b.getHELM()
        print("The HELM format is: {}".format(helm))
        logger.info("Reading biln molecule: {}".format(args.biln))
        logger.info("The helm molecule is: {}".format(helm))

    # 2. A single HELM input directly in the command line
    elif args.helm:
        b = BILN(helm=args.helm)
        biln = b.getBILN()
        print(args.helm)
        print("The BILN format is: {}".format(biln))
        logger.info("Reading helm molecule: {}".format(args.helm))
        logger.info("The biln molecule is: {}".format(biln))

    # 3. A file with a list of BILN molecules
    elif args.table_biln:
        report = open('report_biln.txt','w')
        with open(args.table_biln) as inTable:
            for line in inTable:
                try:
                    b = BILN(biln=line.strip())
                    helm = b.getHELM()
                    report.write(helm+'\n')
                    logger.info("Reading biln molecule: {}".format(line))
                except Exception as e:
                    logger.warning("Failed to process {},{}".format(line,e))
        print("The HELM formats were generated and saved in report_biln.txt")

    # 4. A file with a list of HELM molecules
    elif args.table_helm:
        report = open('report_helm.txt','w')
        with open(args.table_helm) as inTable:
            for line in inTable:
                try:
                    b = BILN(helm=line.strip())
                    biln = b.getBILN()
                    report.write(biln+'\n')
                    logger.info("Reading helm molecule: {}".format(line))
                except Exception as e:
                    logger.warning("Failed to process {},{}".format(line,e))
        print("The BILN formats were generated and saved in report_helm.txt")

    # Termination and cleanup.
    logger.info("Successful completion")
