"""This module supports the creation of AUIDs (area unique identifiers) for DGGSs.
  It is adapted from the code in <https://github.com/IAAA-Lab/dggs-auids>. The implemented algorithms are
  described in:
  "On the problem of providing unique identifiers for areas with any shape on Discrete Global Grid Systems"
  R. BÉJAR, M.Á. LATRE, F.J. LOPEZ-PELLICER, J. NOGUERAS-ISO, F.J. ZARAZAGA-SORIA.
  Accepted Short Papers and Posters from the 22nd AGILE Conference on Geo-information Science.
  Cyprus University of Technology 17-20 June 2019, Limassol, Cyprus, 2019.
  PDF: <https://agile-online.org/images/conferences/2019/documents/short_papers/58_Upload_your_PDF_file.pdf>
"""

import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())