ProFeatMap is a Python based program creating highly customizable 2D maps showing elements of interests (domains, repeats, PTM sites,...) called feature as defined in Uniprot. Several shapes, all hexcolors and more than a hundred colormaps, are available to modify most elements on the map. Numerical values can be displayed on one or several features of interest. ProFeatMap also compiles all features extracted from the Uniprot database in a single file: List of all features, occurence of each feature, length of proteins, sequence of proteins and a list of all available PDB structures. ProFeatMap can also extract sequences for a chosen feature, usable for multiple alignment purposes. It is also possible to search for motifs or other conserved elements using regular expressions.



You can access it here: profeatmap.pythonanywhere.com/

Or you can install it locally by following these steps:

1/ Get a copy of all files from the github repository
2/ Install Python 3.9.5
3/ Install following packages :
- openpyxl 3.0.7
- odfpy 1.4.1
- dash 2.0.0
- dash-bootstrap-components 1.0.3
- dash_extensions 0.0.60
- numpy 1.19.1
- pandas 1.4.1
- urllib3 1.25.9
- pillow 7.2.0
- matplotlib 3.5.1
- xlsxwriter 1.2.9
- uuid 1.30
4/ Run this app with `python app.py` and visit http://127.0.0.1:8050/ in your web browser.

