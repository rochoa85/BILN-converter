# BILN-converter

Script to manipulate BILN format and convert it to HELM (and viceversa)

## Basics

The “Boehringer Ingelheim Line Notation” (BILN) is a human-readable line format that can describe even complex peptides containing multiple chemical modifications. Similar to HELM and other approaches, BILN connects the high-level line notation and the underlying atomic description with a monomer library. A monomer consists of atoms and bonds and can be defined by any chemical structure format (e.g. SMILES or an SD MolBlock). It has a unique identifier used in the BILN and one or several attachment points which define the possible connections to adjacent monomers. Because of that, we allow an easy interconversion between the formats depending on the user needs.

In the following examples, we provide some basic command line codes to convert single or multiple sequences between the two formats. The files used are provided within the code

## Examples

The main script `BILN.py` has the following syntax:

```
usage: BILN.py [-h]
               (--biln text | --helm text | --table_biln filename | --table_helm filename)
               [--logfile filename] [-v]
```

The specific arguments are:

```
  -h, --help            show this help message and exit
  --biln text           BILN string to convert into HELM format.
  --helm text           HELM string to convert into BILN format.
  --table_biln filename
                        Table containing BILN strings to convert into HELM.
  --table_helm filename
                        Table containing HELM strings to convert into BILN.

Logging options:
  --logfile filename    Output messages to given logfile, default is stderr.
  -v, --verbose         Increase output verbosity

```

A first example using a single BILN string can be run as:

`python BILN.py --biln 'Ac(1,2).A-K(1,3)(2,2).Me(2,1)'`

The output will be printed in the shell based on the function logging:

```
INFO:Reading biln molecule: Ac(1,2).A-K(1,3)(2,2).Me(2,1)
INFO:The helm molecule is: PEPTIDE1{[Ac]}|PEPTIDE2{A.K}|PEPTIDE3{[Me]}$PEPTIDE1,PEPTIDE2,1:R2-2:R3|PEPTIDE2,PEPTIDE3,2:R2-1:R1$$$V2.0
INFO:Successful completion
```

The HELM notation can be obtained from the output or sent to a file. Similarly, a HELM string can be used as input to generate the corresponding BILN format. However, if we have a list of HELM strings, it is possible to provide a text file containing all the sequences to be converted (`list_helm.txt`):

```
PEPTIDE1{[Ac]}|PEPTIDE2{A.K}$PEPTIDE1,PEPTIDE2,1:R2-2:R3$$$V2.0
PEPTIDE1{[Ac]}|PEPTIDE2{A.K}|PEPTIDE3{[Me]}$PEPTIDE1,PEPTIDE2,1:R2-2:R3|PEPTIDE2,PEPTIDE3,2:R2-1:R1$$$V2.0
PEPTIDE1{D.T.H.F.P.I.C.I.F.C.C.G.C.C.H.R.S.K.C.G.M.C.C.K.T}$PEPTIDE1,PEPTIDE1,7:R3-23:R3|PEPTIDE1,PEPTIDE1,10:R3-13:R3|PEPTIDE1,PEPTIDE1,11:R3-19:R3|PEPTIDE1,PEPTIDE1,14:R3-22:R3$$$V2.0
PEPTIDE1{F.V.N.Q.H.L.C.G.S.H.L.V.E.A.L.Y.L.V.C.G.E.R.G.F.F.Y.T.P.K.T}|PEPTIDE2{G.I.V.E.Q.C.C.T.S.I.C.S.L.Y.Q.L.E.N.Y.C.N}$PEPTIDE1,PEPTIDE2,7:R3-7:R3|PEPTIDE1,PEPTIDE2,19:R3-20:R3|PEPTIDE2,
PEPTIDE2,6:R3-11:R3$$$V2.0
PEPTIDE1{H.[Aib].E.G.T.F.T.S.D.V.S.S.Y.L.E.G.Q.A.A.K.E.F.I.A.W.L.V.R.G.R.G}|PEPTIDE2{[C18DA].[gGlu].[OEG].[OEG]}$PEPTIDE1,PEPTIDE2,20:R3-4:R2$$$V2.0
PEPTIDE1{[Abu].[Sar].[NMeL].V.[NMeL].A.[DAla].[NMeL].[NMeL].[NMeV].[NMeThr4RBut2enyl]}$PEPTIDE1,PEPTIDE1,1:R1-11:R2$$$V2.0
```

The command in this case can be called as:

`python BILN.py --table_helm list_helm.txt`

The output will be a file called `report_helm.txt` with the BILN conversions:

```
Ac(1,2).A-K(1,3)
Ac(1,2).A-K(1,3)(2,2).Me(2,1)
D-T-H-F-P-I-C(1,3)-I-F-C(2,3)-C(3,3)-G-C(2,3)-C(4,3)-H-R-S-K-C(3,3)-G-M-C(4,3)-C(1,3)-K-T
F-V-N-Q-H-L-C(1,3)-G-S-H-L-V-E-A-L-Y-L-V-C(2,3)-G-E-R-G-F-F-Y-T-P-K-T.G-I-V-E-Q-C(3,3)-C(1,3)-T-S-I-C(3,3)-S-L-Y-Q-L-E-N-Y-C(2,3)-N
H-Aib-E-G-T-F-T-S-D-V-S-S-Y-L-E-G-Q-A-A-K(1,3)-E-F-I-A-W-L-V-R-G-R-G.C18DA-gGlu-OEG-OEG(1,2)
Abu(1,1)-Sar-NMeL-V-NMeL-A-DAla-NMeL-NMeL-NMeV-NMeThr4RBut2enyl(1,2)
```

## Support

For more information about BILN please refer to the publication: 'BILN – A Human-readable Line Notation for Complex Peptides', JCIM, 2022.
