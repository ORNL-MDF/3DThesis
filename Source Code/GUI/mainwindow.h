//This software has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy. 
//Research was co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) – Transformer Resilience and Advanced Components (TRAC) Program.

/*Copyright 2019 UT-Battelle, LLC
*
* All Rights Reserved
*
* Authors: Benjamin Stump <stumpbc@ornl.gov>, Alex Plotkowski, James Ferguson, Kevin Sisco
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice,
*	 this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. Neither the name of 3DThesis nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QStyleFactory>
#include <QDir>
#include <QFileDialog>
#include <QTextStream>
#include <QtConcurrent/QtConcurrent>
#include <string>
#include <string>
#include <fstream>
#include <iostream>
#include <QMessageBox>

#include "modelwrapper.h"
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    ModelWrapper * Model;


public slots:
    void UpdateProgress(int progress);

private slots:

    // Start Button For Simulation
    void on_RunSimulation_clicked();

    // Opens File Dialog for Output Directory
    void on_BrowseOutput_clicked();

    // Opens File Dialog for Path File
    void on_BrowsePath_clicked();

    // Opens File Dialog for Importing Previous Simulation
    void on_Import_triggered();


    // Updating Functions.. These need to be Grouped.
    void on_xmin_editingFinished();
    void on_ymin_editingFinished();
    void on_zmin_editingFinished();
    void on_xmax_editingFinished();
    void on_ymax_editingFinished();
    void on_zmax_editingFinished();
    void on_imax_editingFinished();
    void on_jmax_editingFinished();
    void on_kmax_editingFinished();

    // Stop Button for Simulation
    void on_StopSim_clicked();

private:
    Ui::MainWindow *ui;

    // Function to Generate output folders and SimParams.txt
    void CreateOutputDestination();

    // Loads in Path File and Checks if it is correct/Sets maximum and minimum values
    bool BeamPathFileCheck(QString BeamPathFile);

    // Updates the spatial resoultion of the Simulation
    void UpdateResolution();

    //Checks if all simulation fields are filled out.
    QString CheckingSimulationValues();
};

#endif // MAINWINDOW_H
