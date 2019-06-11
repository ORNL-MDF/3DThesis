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

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    // Setting Style for the GUI
    QApplication::setStyle(QStyleFactory::create("Fusion"));
    ui->setupUi(this);
    Model = new ModelWrapper();

    connect(Model,SIGNAL(progress(int)),this,SLOT(UpdateProgress(int)));
}

MainWindow::~MainWindow()
{
    delete ui;
    delete Model;
}

void MainWindow::on_RunSimulation_clicked()
{
    ui->RunSimulation->setEnabled(false);
    QString Check = CheckingSimulationValues();
    if( Check.length() > 4 )
    {
        QString Message = "A Simulation Property is Blank in " + Check;
        QMessageBox Box;
        Box.setText(Message);
        Box.exec();
        ui->RunSimulation->setEnabled(true);
        return;

    }


    CreateOutputDestination();

    QString Main,Path;

    // Defines SimParams to Read in
    Main = (ui->OutputDir->text() + "/" + ui->SimName->text() + "/Data/" + ui->SimName->text() + ".txt");
    Path = ui->PathDir->text();

    Model->setPaths(Main,Path);
    Model->start();


}

void MainWindow::CreateOutputDestination()
{
    // Creating Directories for Simulation
    qDebug() << "OutputString:" << (ui->OutputDir->text() + ui->SimName->text());
    QDir().mkdir((ui->OutputDir->text() + "/" + ui->SimName->text())); // Creating Main Directory
    QDir().mkdir((ui->OutputDir->text() + "/" + ui->SimName->text() + "/Data/"));
    QDir().mkdir((ui->OutputDir->text() + "/" + ui->SimName->text() + "/Data/" + ui->SimName->text() + "(Data)/"));



    // Generating Simulation File for Running Simulation
    QFile File(ui->OutputDir->text() + "/" + ui->SimName->text() + "/Data/" + ui->SimName->text() + ".txt");
    if (File.open(QIODevice::ReadWrite)) {
        QTextStream stream(&File);
        stream << ui->SimName->text() << "  //Sim Name" << endl;
        stream << ui->OutputDir->text() + "/" + ui->SimName->text() << "  //Path to output folder (Can be relative or absolute) Folder must already exist!!" << endl;
        stream << ui->NumThreads->text() << "  //number of parallel threads (thnum)" << endl;
        stream << ui->imax->text() << " " <<ui->jmax->text() << " " << ui->kmax->text() << "  //imax,jmax,kmax" << endl;
        stream << ui->xmin->text().toDouble()/1000 << " " << ui->xmax->text().toDouble()/1000 << "  // xmin, xmax [m]" << endl;
        stream << ui->ymin->text().toDouble()/1000 << " " << ui->ymax->text().toDouble()/1000 << "  // ymin, ymax [m]" << endl;
        stream << ui->zmin->text().toDouble()/1000 << " " << ui->zmax->text().toDouble()/1000 << "  // zmin, zmax [m]" << endl;
        stream << ui->SimMode->text() << "  //mode" << endl;
        stream << ui->SimTimeStep->text() << " //Time step(dt) [s]" << endl;
        stream << ui->OutputFreq->text() << "  //Out_freq" << endl;
        stream << ui->PreheatTemp->text() << " " << ui->SolidTemp->text() << "Preheat and Solidification Temperature [k]" << endl;
        stream << ui->SpecificHeat->text() << "  //Specific Heat [J/kg-K]" << endl;
        stream << ui->ThermalConductivity->text() << "  //Thermal conductivity [W/mk]" << endl;
        stream << ui->Density->text() << "  //Density [kg/m3]" << endl;
        stream << (ui->BeamAx->text().toDouble())/(1000000) << " " << (ui->BeamAy->text().toDouble())/(1000000) << " " << (ui->BeamAz->text().toDouble())/(1000000) << " //Beam Width(ax,ay,ax) [m]" << endl;
        stream << (ui->AbsorptionEff->text().toDouble()/100) << "  // Absportion efficiency (eff)" << endl;;
        stream << ui->BeamPower->text() << "  //Beam powder(q) [W]" << endl;
    File.close();
    }

}

void MainWindow::on_BrowseOutput_clicked()
{

    // Creating a File Dialog to Check where to put the Simulation Output Folders
    QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),"/home",QFileDialog::ShowDirsOnly| QFileDialog::DontResolveSymlinks);
    if(dir.isNull())
        return;
    ui->OutputDir->setText(dir);
}

void MainWindow::on_BrowsePath_clicked()
{
    // Creating a File Dialog to Check where the Simulation Path is
    QString Pathfilename = QFileDialog::getOpenFileName(this,tr("Open Document"), QDir::currentPath(), tr("text files (*.txt)"));
    if(Pathfilename.isNull())
        return;


    // Checking if Path File is Valid, as well as Finding the Maximum and Minimum Values in the X,Y, and Z directions
    bool Check = false;
    Check = BeamPathFileCheck(Pathfilename);
    if(Check == true)
        ui->PathDir->setText(Pathfilename);


    // Defining the Resolutions of the X,Y,and Z Directions
    UpdateResolution();
}

bool MainWindow::BeamPathFileCheck(QString BeamPathFile)
{


    // Reading In Beam Path File to Determine max and min values for simulation
    QFile File(BeamPathFile);
    QTextStream in(&File);
    QStringList List;
    File.open(QIODevice::ReadOnly);
    do {
        QString line = in.readLine();
        List.append(line);
    } while(!in.atEnd());
    File.close();

    QVector <QVector<double>> Points;
    // Reading in Path File into QStringList
    for(int i = 1; i < List.length(); i++)
    {
         QStringList list = List[i].split(QRegExp("\\s"),QString::SkipEmptyParts);
         // Removing Skipped Lines
         if(list.length()==0)
             continue;
         // Basic Way of Checking Formatting
         if(list.length()!=6)
             return false;

         // Putting List in Vector to Get max and min values.
        for(int j=0;j<list.length();j++)
        {
            QVector<double> Temp;
            Temp.append(list[1].toDouble());
            Temp.append(list[2].toDouble());
            Temp.append(list[3].toDouble());
            Points.append(Temp);
        }
    }

    // Pre defining values to start max/min search
    double xmin, xmax, ymin, ymax, zmin, zmax;
    xmax = Points[0][0]; xmin = Points[0][0];
    ymin = Points[0][1]; ymax = Points[0][1];
    zmin = Points[0][2]; zmax = Points[0][2];


    // Max/Min search
    for(int i = 0; i < Points.length(); i++)
    {
        double x,y,z;
        x = Points[i][0]; y = Points[i][1]; z = Points[i][2];

        if(x>xmax) // Checking if Value is larger than Max
        {xmax = x;}
        else if( x < xmin) // Checking if Value is smaller than min
        {xmin =x;}

        if(y>ymax) // Checking if Value is larger than Max
            ymax = y;
        else if( y < ymin) // Checking if Value is smaller than min
        {ymin =y;}

        if(z>zmax) // Checking if Value is larger than Max
            zmax = z;
        else if( z < zmin) // Checking if Value is smaller than min
        {zmin =z;}
    }



    // Setting GUI QlineEdits to max and min values
    ui->xmin->setText((QString::number(xmin - 1)));
    ui->xmax->setText((QString::number(xmax + 1)));
    ui->ymin->setText((QString::number(ymin - 1)));
    ui->ymax->setText((QString::number(ymax + 1)));
    ui->zmin->setText((QString::number(zmin - .05)));
    ui->zmax->setText((QString::number(zmax)));
    return true;


}

void MainWindow::UpdateResolution()
{
    // This Function Calculates the Spatial resolution of the Simulation Based on Input Paramaters.
    double i, j ,k;
    double xmax, xmin;
    double ymax, ymin;
    double zmax, zmin;
    i    = ui->imax->text().toDouble(); j    = ui->jmax->text().toDouble();  k    = ui->kmax->text().toDouble();
    xmax = ui->xmax->text().toDouble(); xmin = ui->xmin->text().toDouble();  ymax = ui->ymax->text().toDouble();
    ymin = ui->ymin->text().toDouble(); zmax = ui->zmax->text().toDouble();  zmin = ui->zmin->text().toDouble();
    ui->xRes->setText(QString::number(std::abs(xmax-xmin)/(i-1)*1000));
    ui->yRes->setText(QString::number(std::abs(ymax-ymin)/(j-1)*1000));
    ui->zRes->setText(QString::number(std::abs(zmax-zmin)/(k-1)*1000));

}

void MainWindow::on_xmin_editingFinished()
{
    // Checkign for When Xmin Changes to Update Resolution of Simulation
    UpdateResolution();
}

void MainWindow::on_ymin_editingFinished()
{
    // Checkign for When ymin Changes to Update Resolution of Simulation
    UpdateResolution();
}

void MainWindow::on_zmin_editingFinished()
{
    // Checkign for When zmin Changes to Update Resolution of Simulation
    UpdateResolution();
}

void MainWindow::on_xmax_editingFinished()
{
    // Checkign for When Xmax Changes to Update Resolution of Simulation
    UpdateResolution();
}

void MainWindow::on_ymax_editingFinished()
{
    // Checkign for When ymax Changes to Update Resolution of Simulation
    UpdateResolution();
}

void MainWindow::on_zmax_editingFinished()
{
    // Checkign for When zmax Changes to Update Resolution of Simulation
    UpdateResolution();
}

void MainWindow::on_imax_editingFinished()
{
    // Checkign for When imax Changes to Update Resolution of Simulation
    UpdateResolution();
}

void MainWindow::on_jmax_editingFinished()
{
    // Checkign for When jmax Changes to Update Resolution of Simulation
    UpdateResolution();
}

void MainWindow::on_kmax_editingFinished()
{
    // Checkign for When kmax Changes to Update Resolution of Simulation
    UpdateResolution();
}

void MainWindow::on_Import_triggered()
{

    QString Pathfilename = QFileDialog::getOpenFileName(this,tr("Open Document"), QDir::currentPath(), tr("text files (*.txt)"));
    if(Pathfilename.isNull())
        return;

    std::string line;
    std::ifstream simfile;
    simfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);


    //Defining Variables for Read input
    std::string sim_name; //Name for output data files
    std::string out_path; //Directory of output data files
    //Material properties
    double kon, rho, cps, Tsol, Tinit;
    //Beam parameters
    double ax, ay, az;
    double eff, q;
    //Simulation parameters
    int thnum;
    int imax, jmax, kmax;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    int mode, out_freq;
    double dt;

    std::string path = Pathfilename.toStdString();
    int count = 0;
    try{
        simfile.open(path, std::ios::in);
        simfile >> sim_name; std::getline(simfile, line); count++;
        simfile >> out_path; std::getline(simfile, line); count++;
        simfile >> thnum; getline(simfile, line); count++;	 		//Number of threads to use in parallel
        simfile >> imax >> jmax >> kmax; getline(simfile, line); count++;	//Number of points in each direction
        simfile >> xmin >> xmax; getline(simfile, line); count++;	//x-direction limits
        simfile >> ymin >> ymax; getline(simfile, line); count++;   //y-direction limits
        simfile >> zmin >> zmax; getline(simfile, line); count++;   //z-direction limits
        simfile >> mode; getline(simfile, line); count++;   //Simulation mode (1 - normal, 2 - GV only)
        simfile >> dt; getline(simfile, line); count++;	//Time step
        simfile >> out_freq; getline(simfile, line); count++;	//Output frequency in number of time steps
        //material part
        simfile >> Tinit >> Tsol; getline(simfile, line); count++; //Pre-heat and Solidification Temperatures
        simfile >> cps; getline(simfile, line);	count++;//Specific heat
        simfile >> kon; getline(simfile, line);	count++;//Thermal conductivity
        simfile >> rho; getline(simfile, line);	count++;//Density
        //beam size
        simfile >> ax >> ay >> az; getline(simfile, line); count++;//Beam widths
        simfile >> eff; getline(simfile, line); count++;		//Efficiency of absorption
        simfile >> q; getline(simfile, line); count++;		//Power

        simfile.close();
    }
    catch (const std::ifstream::failure& e){
        std::cout << "Exception opening/reading sim file at " << count << std::endl;
    }

    // Setting GUI LineEdits to the Read in Values
    ui->SimName->setText(QString::fromStdString(sim_name));
    ui->NumThreads->setText(QString::number(thnum));
    ui->imax->setText(QString::number(imax)); ui->jmax->setText(QString::number(jmax)); ui->kmax->setText(QString::number(kmax));
    ui->xmin->setText(QString::number((xmin*1000))); ui->ymin->setText(QString::number((ymin*1000))); ui->zmin->setText(QString::number((zmin*1000)));
    ui->xmax->setText(QString::number((xmax*1000))); ui->ymax->setText(QString::number((ymax*1000))); ui->zmax->setText(QString::number((zmax*1000)));
    ui->SimMode->setText(QString::number(mode)); ui->SimTimeStep->setText(QString::number(dt)); ui->OutputFreq->setText(QString::number(out_freq));
    ui->SolidTemp->setText(QString::number(Tsol)); ui->PreheatTemp->setText(QString::number(Tinit));
    ui->SpecificHeat->setText(QString::number(cps)); ui->ThermalConductivity->setText(QString::number(kon)); ui->Density->setText(QString::number(rho));
    ui->BeamAx->setText(QString::number(ax*1000000)); ui->BeamAy->setText(QString::number(ay*1000000)); ui->BeamAz->setText(QString::number(az*1000000));
    ui->BeamPower->setText(QString::number(q)); ui->AbsorptionEff->setText(QString::number(eff*100));
    UpdateResolution();

}

void MainWindow::on_StopSim_clicked()
{
    // Killing Simulation Thread on Stop button clicked.. This is a non-safe thread exit and needs to be redone.
    Model->terminate();
    ui->RunSimulation->setEnabled(true);
}

void MainWindow::UpdateProgress(int progress)
{
    ui->progressBar->setValue(progress);

    if(progress == 100)
        ui->RunSimulation->setEnabled(true);
}

QString MainWindow::CheckingSimulationValues()
{
    // File Inputs Check
    if(ui->OutputDir->text().isEmpty() | ui->PathDir->text().isEmpty() | ui->SimName->text().isEmpty())
        return "File Inputs";
    // Material Properties Check
    if(ui->PreheatTemp->text().isEmpty() | ui->SolidTemp->text().isEmpty() | ui->SpecificHeat->text().isEmpty() | ui->ThermalConductivity->text().isEmpty() | ui->Density->text().isEmpty())
        return "Material Prop.";
    // Laser Properties Check
    if(ui->BeamAx->text().isEmpty() | ui->BeamAy->text().isEmpty() | ui->BeamAz->text().isEmpty() | ui->AbsorptionEff->text().isEmpty() | ui->BeamPower->text().isEmpty())
        return "Laser Prop.";
    // Simulation Bounds Check
    if(ui->xmin->text().isEmpty() | ui->ymin->text().isEmpty() | ui->zmin->text().isEmpty() | ui->xmax->text().isEmpty() | ui->ymax->text().isEmpty() |  ui->zmax->text().isEmpty() | ui->imax->text().isEmpty() | ui->jmax->text().isEmpty() | ui->kmax->text().isEmpty())
        return "Simulation Bounds";
    // Settings Check
    if(ui->NumThreads->text().isEmpty() | ui->SimMode->text().isEmpty() | ui->OutputFreq->text().isEmpty())
        return "Settings";


    return "Good";


}
