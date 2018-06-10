# -*- coding: utf-8 -*-
"""
Created on Sun Jun 10 17:30:49 2018

@author: anton
"""  
    
import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QVBoxLayout, QPushButton, QLineEdit, QMessageBox
import GEOparse
import numpy as np
import pandas as pd
from scipy.stats.stats import ks_2samp
 
def GSEA (geo_ID, gene_list):
    gse = GEOparse.get_GEO(geo=geo_ID, destdir="./")
    expression = gse.pivot_samples('VALUE').T
    experiments = {}
    for i, (idx, row) in enumerate(gse.phenotype_data.iterrows()):
        tmp = {}
        tmp["Type"] = 1 if "control" in row["description"] else 0
        experiments[i] = tmp
    experiments = pd.DataFrame(experiments).T
    counter = 0
    all_genes_set = []
    all_corr_set = []
    genes_corr_set = []
    for gene in expression:
        counter += 1
        if counter <= 3:
            continue
        all_genes_set.append(gene)               
        corr_matrix = np.corrcoef([list(experiments['Type']), list(expression[gene])])
        all_corr_set.append(corr_matrix[0,1])
        if gene in gene_list:
            genes_corr_set.append(corr_matrix[0,1])
    p_value = ks_2samp(genes_corr_set, all_corr_set)[1]
    return(str(p_value))

class App(QMainWindow):
 
    def __init__(self):
        super().__init__()
        self.title = 'GSEA'
        self.left = 10
        self.top = 10
        self.width = 400
        self.height = 200
        self.initUI()
 
    def initUI(self):
        self.layout = QVBoxLayout()
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
 

        self.textbox = QLineEdit(self)
        self.textbox.move(20, 20)
        self.textbox.resize(280,40)
        self.textbox.setToolTip("Please enter GEO id")
        self.textbox.resize(self.textbox.sizeHint())
        
        self.textbox2 = QLineEdit(self)
        self.textbox2.move(20, 90)
        self.textbox2.resize(280,40)
        self.textbox2.setToolTip("Please enter Gene list")
        self.textbox2.resize(self.textbox2.sizeHint())


        self.button = QPushButton('Get p_value', self)
        self.button.move(20,150)
 
        self.button.clicked.connect(self.on_click)
        self.show()
    

    def on_click(self):
        self.geo_ID = self.textbox.text()
        self.gene_list = self.textbox2.text().split()
        QMessageBox.question(self, 'GSEA p_value', "p_value: " + GSEA(self.geo_ID,self.gene_list), QMessageBox.Ok, QMessageBox.Ok)
        self.textbox.setText("")
 
 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    app.exec_()
