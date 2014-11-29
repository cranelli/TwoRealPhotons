""" HTML Table - creates an HTML table from added cells """

# Python Module For making an HTML Table
# Created by Nicholas Zube, 8.22.2014
# Tested in Python 2.6.4

class HTMLTable:
    
    def __init__(self, rows = 1, cols = 1, fileName = "HTMLTable"):       
        self.nRows = 1
        self.nCols = 1
        self.fileTitle = fileName       
        if(rows > 0): self.nRows = rows;
        if(cols > 0): self.nCols = cols;
        
        self.cells = []
        self.cellAlign = []    
        self.cellsSingleRow = [""] * self.nCols
        self.cellAlignSingleRow = ["middle"] * self.nCols
        for i in xrange(self.nRows):
            self.cells.append(self.cellsSingleRow[:])
            self.cellAlign.append(self.cellAlignSingleRow[:])
        self.rowAlign = ["center"] * self.nRows
        self.tSettings = ["title", "2", "2", "3", "border", "all", 
                          "summary of table"]
        self.currentRow = 0
        self.currentCol = 0


    def BuildFromFile(self, fileLocation):
        with open(fileLocation) as infile:
            HTMLFileName = fileLocation.translate(None,".txt") + ".html"          
            bRows = 0
            bCols = 0
            bcells = []
            delimiter = "\t"
            for line in infile:
                # Create a list of "cell strings" separated by the delimiter
                lineList = line.split(delimiter)
                # Remove \n from line string
                for entry in xrange(len(lineList)):
                    lineList[entry] = lineList[entry].translate(None, "\n")
                if len(lineList) > bCols: 
                    bCols = len(lineList)
                    
                bcells.append(lineList)
                bRows += 1                               
            table = HTMLTable(bRows, bCols)            
            
            # Add blank strings to make all rows same size as largest
            for row in bcells:
                while len(row) < bCols:
                    row.append("-")
            # Fill table
            for row in xrange(bRows):
                for col in xrange(bCols):
                    table.Add(bcells[row][col], row, col)
                                        
            table.Print(HTMLFileName)         

    def Add(self, entry = "", row = 0, col = 0):
        self.GoTo(row,col)
        self.cells[self.currentRow][self.currentCol] = entry
        # Move the next horizontal cell
        self.currentCol += 1        
        # If not on last row
        if self.currentRow < self.nRows-1:
            # If beyond final column in row, go to beginning of next row
            if self.currentCol >= self.nCols:
                self.currentRow += 1
                self.currentCol = 0
        # Else is on last row; if beyond final column, return to final    
        #elif self.currentCol >= self.nCols:
            #self.nCols -= 1

    def GetCurrentRow(self):
        return(self.currentRow)

    def GetCurrentCol(self):
        return(self.currentCol)

    def GoTo(self, row, col):
        if abs(row) < self.nRows: 
            self.currentRow = abs(row)
        # Else beyond max, so set to max
        else: self.currentRow = self.nRows-1
        if abs(col) < self.nCols: 
            self.currentCol = abs(col)
        else: self.currentCol = self.nCols-1
        
    def Look(self, row=-1, col=-1):
        if row<self.nRows and col<self.nCols and row>=0 and col>=0:
            co = (self.cells[row][col] + ", " + str(row) + ", " + str(col))
            return co
        else: 
            co = (self.cells[self.GetCurrentRow()][self.GetCurrentCol()-1]
                 + ", " + str(self.GetCurrentRow()) + ", "
                 + str(self.GetCurrentCol()-1))
            return (co, "(", str(self.nRows), str(self.nCols),")")
        
    def Output(self):
        output = ("<table id=\"{0}\" border=\"{1}\" cellspacing=\"{2}\" cellpadding=\"{3}\" frame=\"{4}\" rules=\"{5}\" summary=\"{6}\"><tbody>"
                  ).format(*tuple(self.tSettings))
        for row in range(self.nRows):
            output += "</tr><tr align=\"{0}\">".format(self.rowAlign[row])
            for col in range(self.nCols):
                output += "<td valign=\"{0}\">{1}</td>".format(self.cellAlign[row][col], self.cells[row][col])
        output += "</tr></tbody></table>\n"
        return output

    def Print(self, outfileName = "HTMLTable.html"):
        with open(outfileName,'w') as outfile:          
            outfile.write(self.Output())

    def SetParams(self, tid = "title", border = "2", cellspacing = "2", 
                  cellpadding = "3", frame = "border", rules = "all", 
                  summary = "summary of table"):
        self.tSettings = [tid, border, cellspacing, cellpadding, frame, 
                          rules, summary]
