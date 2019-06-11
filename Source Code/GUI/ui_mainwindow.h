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

/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.10.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *Import;
    QWidget *centralWidget;
    QVBoxLayout *verticalLayout_10;
    QHBoxLayout *horizontalLayout;
    QVBoxLayout *verticalLayout;
    QLabel *SimNameLabel;
    QLabel *OutputLabel;
    QLabel *PathLabel;
    QVBoxLayout *verticalLayout_2;
    QLineEdit *SimName;
    QLineEdit *OutputDir;
    QLineEdit *PathDir;
    QVBoxLayout *verticalLayout_3;
    QSpacerItem *verticalSpacer_5;
    QSpacerItem *horizontalSpacer;
    QPushButton *BrowseOutput;
    QPushButton *BrowsePath;
    QFrame *line;
    QLabel *label_37;
    QTabWidget *tabWidget;
    QWidget *tab;
    QVBoxLayout *verticalLayout_15;
    QHBoxLayout *horizontalLayout_6;
    QVBoxLayout *verticalLayout_4;
    QLabel *label_6;
    QLabel *label_7;
    QLabel *label_8;
    QLabel *label_9;
    QVBoxLayout *verticalLayout_5;
    QLineEdit *SolidTemp;
    QLineEdit *SpecificHeat;
    QLineEdit *ThermalConductivity;
    QLineEdit *Density;
    QVBoxLayout *verticalLayout_6;
    QLabel *label_11;
    QLabel *label_12;
    QLabel *label_13;
    QLabel *label_14;
    QSpacerItem *verticalSpacer_3;
    QWidget *tab_3;
    QVBoxLayout *verticalLayout_11;
    QHBoxLayout *horizontalLayout_7;
    QVBoxLayout *verticalLayout_7;
    QLabel *label_15;
    QLabel *label_16;
    QLabel *label_17;
    QLabel *label_18;
    QLabel *label_19;
    QLabel *label_5;
    QVBoxLayout *verticalLayout_8;
    QLineEdit *BeamAx;
    QLineEdit *BeamAy;
    QLineEdit *BeamAz;
    QLineEdit *AbsorptionEff;
    QLineEdit *BeamPower;
    QLineEdit *PreheatTemp;
    QVBoxLayout *verticalLayout_9;
    QLabel *label_20;
    QLabel *label_21;
    QLabel *label_22;
    QLabel *label_23;
    QLabel *label_24;
    QLabel *label_10;
    QSpacerItem *verticalSpacer_2;
    QWidget *tab_4;
    QVBoxLayout *verticalLayout_30;
    QHBoxLayout *horizontalLayout_9;
    QHBoxLayout *horizontalLayout_3;
    QVBoxLayout *verticalLayout_23;
    QVBoxLayout *verticalLayout_14;
    QLabel *label_25;
    QLabel *label_26;
    QLabel *label_27;
    QVBoxLayout *verticalLayout_19;
    QLabel *label_31;
    QLabel *label_32;
    QLabel *label_33;
    QVBoxLayout *verticalLayout_24;
    QVBoxLayout *verticalLayout_17;
    QLineEdit *xmin;
    QLineEdit *ymin;
    QLineEdit *zmin;
    QVBoxLayout *verticalLayout_20;
    QLineEdit *imax;
    QLineEdit *jmax;
    QLineEdit *kmax;
    QVBoxLayout *verticalLayout_28;
    QLabel *label_44;
    QLabel *label_45;
    QLabel *label_46;
    QLabel *label_47;
    QLabel *label_48;
    QLabel *label_49;
    QFrame *line_2;
    QVBoxLayout *verticalLayout_29;
    QHBoxLayout *horizontalLayout_4;
    QVBoxLayout *verticalLayout_16;
    QLabel *label_28;
    QLabel *label_29;
    QLabel *label_30;
    QVBoxLayout *verticalLayout_18;
    QLineEdit *xmax;
    QLineEdit *ymax;
    QLineEdit *zmax;
    QVBoxLayout *verticalLayout_26;
    QLabel *label_38;
    QLabel *label_39;
    QLabel *label_43;
    QHBoxLayout *horizontalLayout_5;
    QVBoxLayout *verticalLayout_21;
    QLabel *label_34;
    QLabel *label_35;
    QLabel *label_36;
    QVBoxLayout *verticalLayout_22;
    QLabel *xRes;
    QLabel *yRes;
    QLabel *zRes;
    QVBoxLayout *verticalLayout_25;
    QLabel *label_40;
    QLabel *label_41;
    QLabel *label_42;
    QSpacerItem *verticalSpacer_4;
    QWidget *tab_2;
    QVBoxLayout *verticalLayout_27;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout_13;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_4;
    QVBoxLayout *verticalLayout_12;
    QLineEdit *NumThreads;
    QLineEdit *SimMode;
    QLineEdit *SimTimeStep;
    QLineEdit *OutputFreq;
    QSpacerItem *verticalSpacer;
    QProgressBar *progressBar;
    QHBoxLayout *horizontalLayout_8;
    QSpacerItem *horizontalSpacer_2;
    QPushButton *RunSimulation;
    QPushButton *StopSim;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(369, 459);
        Import = new QAction(MainWindow);
        Import->setObjectName(QStringLiteral("Import"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        verticalLayout_10 = new QVBoxLayout(centralWidget);
        verticalLayout_10->setSpacing(6);
        verticalLayout_10->setContentsMargins(11, 11, 11, 11);
        verticalLayout_10->setObjectName(QStringLiteral("verticalLayout_10"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        SimNameLabel = new QLabel(centralWidget);
        SimNameLabel->setObjectName(QStringLiteral("SimNameLabel"));

        verticalLayout->addWidget(SimNameLabel);

        OutputLabel = new QLabel(centralWidget);
        OutputLabel->setObjectName(QStringLiteral("OutputLabel"));

        verticalLayout->addWidget(OutputLabel);

        PathLabel = new QLabel(centralWidget);
        PathLabel->setObjectName(QStringLiteral("PathLabel"));

        verticalLayout->addWidget(PathLabel);


        horizontalLayout->addLayout(verticalLayout);

        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        SimName = new QLineEdit(centralWidget);
        SimName->setObjectName(QStringLiteral("SimName"));

        verticalLayout_2->addWidget(SimName);

        OutputDir = new QLineEdit(centralWidget);
        OutputDir->setObjectName(QStringLiteral("OutputDir"));

        verticalLayout_2->addWidget(OutputDir);

        PathDir = new QLineEdit(centralWidget);
        PathDir->setObjectName(QStringLiteral("PathDir"));

        verticalLayout_2->addWidget(PathDir);


        horizontalLayout->addLayout(verticalLayout_2);

        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        verticalSpacer_5 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_3->addItem(verticalSpacer_5);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Preferred, QSizePolicy::Minimum);

        verticalLayout_3->addItem(horizontalSpacer);

        BrowseOutput = new QPushButton(centralWidget);
        BrowseOutput->setObjectName(QStringLiteral("BrowseOutput"));

        verticalLayout_3->addWidget(BrowseOutput);

        BrowsePath = new QPushButton(centralWidget);
        BrowsePath->setObjectName(QStringLiteral("BrowsePath"));

        verticalLayout_3->addWidget(BrowsePath);


        horizontalLayout->addLayout(verticalLayout_3);


        verticalLayout_10->addLayout(horizontalLayout);

        line = new QFrame(centralWidget);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout_10->addWidget(line);

        label_37 = new QLabel(centralWidget);
        label_37->setObjectName(QStringLiteral("label_37"));

        verticalLayout_10->addWidget(label_37);

        tabWidget = new QTabWidget(centralWidget);
        tabWidget->setObjectName(QStringLiteral("tabWidget"));
        tabWidget->setStyleSheet(QStringLiteral("background-color: rgb(240, 240, 240);"));
        tab = new QWidget();
        tab->setObjectName(QStringLiteral("tab"));
        verticalLayout_15 = new QVBoxLayout(tab);
        verticalLayout_15->setSpacing(6);
        verticalLayout_15->setContentsMargins(11, 11, 11, 11);
        verticalLayout_15->setObjectName(QStringLiteral("verticalLayout_15"));
        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        label_6 = new QLabel(tab);
        label_6->setObjectName(QStringLiteral("label_6"));

        verticalLayout_4->addWidget(label_6);

        label_7 = new QLabel(tab);
        label_7->setObjectName(QStringLiteral("label_7"));

        verticalLayout_4->addWidget(label_7);

        label_8 = new QLabel(tab);
        label_8->setObjectName(QStringLiteral("label_8"));

        verticalLayout_4->addWidget(label_8);

        label_9 = new QLabel(tab);
        label_9->setObjectName(QStringLiteral("label_9"));

        verticalLayout_4->addWidget(label_9);


        horizontalLayout_6->addLayout(verticalLayout_4);

        verticalLayout_5 = new QVBoxLayout();
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        SolidTemp = new QLineEdit(tab);
        SolidTemp->setObjectName(QStringLiteral("SolidTemp"));
        SolidTemp->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        SolidTemp->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_5->addWidget(SolidTemp);

        SpecificHeat = new QLineEdit(tab);
        SpecificHeat->setObjectName(QStringLiteral("SpecificHeat"));
        SpecificHeat->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        SpecificHeat->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_5->addWidget(SpecificHeat);

        ThermalConductivity = new QLineEdit(tab);
        ThermalConductivity->setObjectName(QStringLiteral("ThermalConductivity"));
        ThermalConductivity->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        ThermalConductivity->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_5->addWidget(ThermalConductivity);

        Density = new QLineEdit(tab);
        Density->setObjectName(QStringLiteral("Density"));
        Density->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        Density->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_5->addWidget(Density);


        horizontalLayout_6->addLayout(verticalLayout_5);

        verticalLayout_6 = new QVBoxLayout();
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setObjectName(QStringLiteral("verticalLayout_6"));
        label_11 = new QLabel(tab);
        label_11->setObjectName(QStringLiteral("label_11"));

        verticalLayout_6->addWidget(label_11);

        label_12 = new QLabel(tab);
        label_12->setObjectName(QStringLiteral("label_12"));

        verticalLayout_6->addWidget(label_12);

        label_13 = new QLabel(tab);
        label_13->setObjectName(QStringLiteral("label_13"));

        verticalLayout_6->addWidget(label_13);

        label_14 = new QLabel(tab);
        label_14->setObjectName(QStringLiteral("label_14"));

        verticalLayout_6->addWidget(label_14);


        horizontalLayout_6->addLayout(verticalLayout_6);


        verticalLayout_15->addLayout(horizontalLayout_6);

        verticalSpacer_3 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_15->addItem(verticalSpacer_3);

        tabWidget->addTab(tab, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QStringLiteral("tab_3"));
        verticalLayout_11 = new QVBoxLayout(tab_3);
        verticalLayout_11->setSpacing(6);
        verticalLayout_11->setContentsMargins(11, 11, 11, 11);
        verticalLayout_11->setObjectName(QStringLiteral("verticalLayout_11"));
        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        verticalLayout_7 = new QVBoxLayout();
        verticalLayout_7->setSpacing(6);
        verticalLayout_7->setObjectName(QStringLiteral("verticalLayout_7"));
        label_15 = new QLabel(tab_3);
        label_15->setObjectName(QStringLiteral("label_15"));

        verticalLayout_7->addWidget(label_15);

        label_16 = new QLabel(tab_3);
        label_16->setObjectName(QStringLiteral("label_16"));

        verticalLayout_7->addWidget(label_16);

        label_17 = new QLabel(tab_3);
        label_17->setObjectName(QStringLiteral("label_17"));

        verticalLayout_7->addWidget(label_17);

        label_18 = new QLabel(tab_3);
        label_18->setObjectName(QStringLiteral("label_18"));

        verticalLayout_7->addWidget(label_18);

        label_19 = new QLabel(tab_3);
        label_19->setObjectName(QStringLiteral("label_19"));

        verticalLayout_7->addWidget(label_19);

        label_5 = new QLabel(tab_3);
        label_5->setObjectName(QStringLiteral("label_5"));

        verticalLayout_7->addWidget(label_5);


        horizontalLayout_7->addLayout(verticalLayout_7);

        verticalLayout_8 = new QVBoxLayout();
        verticalLayout_8->setSpacing(6);
        verticalLayout_8->setObjectName(QStringLiteral("verticalLayout_8"));
        BeamAx = new QLineEdit(tab_3);
        BeamAx->setObjectName(QStringLiteral("BeamAx"));
        BeamAx->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        BeamAx->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_8->addWidget(BeamAx);

        BeamAy = new QLineEdit(tab_3);
        BeamAy->setObjectName(QStringLiteral("BeamAy"));
        BeamAy->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        BeamAy->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_8->addWidget(BeamAy);

        BeamAz = new QLineEdit(tab_3);
        BeamAz->setObjectName(QStringLiteral("BeamAz"));
        BeamAz->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        BeamAz->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_8->addWidget(BeamAz);

        AbsorptionEff = new QLineEdit(tab_3);
        AbsorptionEff->setObjectName(QStringLiteral("AbsorptionEff"));
        AbsorptionEff->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        AbsorptionEff->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_8->addWidget(AbsorptionEff);

        BeamPower = new QLineEdit(tab_3);
        BeamPower->setObjectName(QStringLiteral("BeamPower"));
        BeamPower->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        BeamPower->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_8->addWidget(BeamPower);

        PreheatTemp = new QLineEdit(tab_3);
        PreheatTemp->setObjectName(QStringLiteral("PreheatTemp"));
        PreheatTemp->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        PreheatTemp->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_8->addWidget(PreheatTemp);


        horizontalLayout_7->addLayout(verticalLayout_8);

        verticalLayout_9 = new QVBoxLayout();
        verticalLayout_9->setSpacing(6);
        verticalLayout_9->setObjectName(QStringLiteral("verticalLayout_9"));
        label_20 = new QLabel(tab_3);
        label_20->setObjectName(QStringLiteral("label_20"));

        verticalLayout_9->addWidget(label_20);

        label_21 = new QLabel(tab_3);
        label_21->setObjectName(QStringLiteral("label_21"));

        verticalLayout_9->addWidget(label_21);

        label_22 = new QLabel(tab_3);
        label_22->setObjectName(QStringLiteral("label_22"));

        verticalLayout_9->addWidget(label_22);

        label_23 = new QLabel(tab_3);
        label_23->setObjectName(QStringLiteral("label_23"));

        verticalLayout_9->addWidget(label_23);

        label_24 = new QLabel(tab_3);
        label_24->setObjectName(QStringLiteral("label_24"));

        verticalLayout_9->addWidget(label_24);

        label_10 = new QLabel(tab_3);
        label_10->setObjectName(QStringLiteral("label_10"));

        verticalLayout_9->addWidget(label_10);


        horizontalLayout_7->addLayout(verticalLayout_9);


        verticalLayout_11->addLayout(horizontalLayout_7);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_11->addItem(verticalSpacer_2);

        tabWidget->addTab(tab_3, QString());
        tab_4 = new QWidget();
        tab_4->setObjectName(QStringLiteral("tab_4"));
        verticalLayout_30 = new QVBoxLayout(tab_4);
        verticalLayout_30->setSpacing(6);
        verticalLayout_30->setContentsMargins(11, 11, 11, 11);
        verticalLayout_30->setObjectName(QStringLiteral("verticalLayout_30"));
        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        verticalLayout_23 = new QVBoxLayout();
        verticalLayout_23->setSpacing(6);
        verticalLayout_23->setObjectName(QStringLiteral("verticalLayout_23"));
        verticalLayout_14 = new QVBoxLayout();
        verticalLayout_14->setSpacing(6);
        verticalLayout_14->setObjectName(QStringLiteral("verticalLayout_14"));
        label_25 = new QLabel(tab_4);
        label_25->setObjectName(QStringLiteral("label_25"));

        verticalLayout_14->addWidget(label_25);

        label_26 = new QLabel(tab_4);
        label_26->setObjectName(QStringLiteral("label_26"));

        verticalLayout_14->addWidget(label_26);

        label_27 = new QLabel(tab_4);
        label_27->setObjectName(QStringLiteral("label_27"));

        verticalLayout_14->addWidget(label_27);


        verticalLayout_23->addLayout(verticalLayout_14);

        verticalLayout_19 = new QVBoxLayout();
        verticalLayout_19->setSpacing(6);
        verticalLayout_19->setObjectName(QStringLiteral("verticalLayout_19"));
        label_31 = new QLabel(tab_4);
        label_31->setObjectName(QStringLiteral("label_31"));

        verticalLayout_19->addWidget(label_31);

        label_32 = new QLabel(tab_4);
        label_32->setObjectName(QStringLiteral("label_32"));

        verticalLayout_19->addWidget(label_32);

        label_33 = new QLabel(tab_4);
        label_33->setObjectName(QStringLiteral("label_33"));

        verticalLayout_19->addWidget(label_33);


        verticalLayout_23->addLayout(verticalLayout_19);


        horizontalLayout_3->addLayout(verticalLayout_23);

        verticalLayout_24 = new QVBoxLayout();
        verticalLayout_24->setSpacing(6);
        verticalLayout_24->setObjectName(QStringLiteral("verticalLayout_24"));
        verticalLayout_17 = new QVBoxLayout();
        verticalLayout_17->setSpacing(6);
        verticalLayout_17->setObjectName(QStringLiteral("verticalLayout_17"));
        xmin = new QLineEdit(tab_4);
        xmin->setObjectName(QStringLiteral("xmin"));
        xmin->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        xmin->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_17->addWidget(xmin);

        ymin = new QLineEdit(tab_4);
        ymin->setObjectName(QStringLiteral("ymin"));
        ymin->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        ymin->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_17->addWidget(ymin);

        zmin = new QLineEdit(tab_4);
        zmin->setObjectName(QStringLiteral("zmin"));
        zmin->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        zmin->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_17->addWidget(zmin);


        verticalLayout_24->addLayout(verticalLayout_17);

        verticalLayout_20 = new QVBoxLayout();
        verticalLayout_20->setSpacing(6);
        verticalLayout_20->setObjectName(QStringLiteral("verticalLayout_20"));
        imax = new QLineEdit(tab_4);
        imax->setObjectName(QStringLiteral("imax"));
        imax->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        imax->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_20->addWidget(imax);

        jmax = new QLineEdit(tab_4);
        jmax->setObjectName(QStringLiteral("jmax"));
        jmax->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        jmax->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_20->addWidget(jmax);

        kmax = new QLineEdit(tab_4);
        kmax->setObjectName(QStringLiteral("kmax"));
        kmax->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        kmax->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_20->addWidget(kmax);


        verticalLayout_24->addLayout(verticalLayout_20);


        horizontalLayout_3->addLayout(verticalLayout_24);

        verticalLayout_28 = new QVBoxLayout();
        verticalLayout_28->setSpacing(6);
        verticalLayout_28->setObjectName(QStringLiteral("verticalLayout_28"));
        label_44 = new QLabel(tab_4);
        label_44->setObjectName(QStringLiteral("label_44"));

        verticalLayout_28->addWidget(label_44);

        label_45 = new QLabel(tab_4);
        label_45->setObjectName(QStringLiteral("label_45"));

        verticalLayout_28->addWidget(label_45);

        label_46 = new QLabel(tab_4);
        label_46->setObjectName(QStringLiteral("label_46"));

        verticalLayout_28->addWidget(label_46);

        label_47 = new QLabel(tab_4);
        label_47->setObjectName(QStringLiteral("label_47"));

        verticalLayout_28->addWidget(label_47);

        label_48 = new QLabel(tab_4);
        label_48->setObjectName(QStringLiteral("label_48"));

        verticalLayout_28->addWidget(label_48);

        label_49 = new QLabel(tab_4);
        label_49->setObjectName(QStringLiteral("label_49"));

        verticalLayout_28->addWidget(label_49);


        horizontalLayout_3->addLayout(verticalLayout_28);


        horizontalLayout_9->addLayout(horizontalLayout_3);

        line_2 = new QFrame(tab_4);
        line_2->setObjectName(QStringLiteral("line_2"));
        line_2->setFrameShape(QFrame::VLine);
        line_2->setFrameShadow(QFrame::Sunken);

        horizontalLayout_9->addWidget(line_2);

        verticalLayout_29 = new QVBoxLayout();
        verticalLayout_29->setSpacing(6);
        verticalLayout_29->setObjectName(QStringLiteral("verticalLayout_29"));
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        verticalLayout_16 = new QVBoxLayout();
        verticalLayout_16->setSpacing(6);
        verticalLayout_16->setObjectName(QStringLiteral("verticalLayout_16"));
        label_28 = new QLabel(tab_4);
        label_28->setObjectName(QStringLiteral("label_28"));

        verticalLayout_16->addWidget(label_28);

        label_29 = new QLabel(tab_4);
        label_29->setObjectName(QStringLiteral("label_29"));

        verticalLayout_16->addWidget(label_29);

        label_30 = new QLabel(tab_4);
        label_30->setObjectName(QStringLiteral("label_30"));

        verticalLayout_16->addWidget(label_30);


        horizontalLayout_4->addLayout(verticalLayout_16);

        verticalLayout_18 = new QVBoxLayout();
        verticalLayout_18->setSpacing(6);
        verticalLayout_18->setObjectName(QStringLiteral("verticalLayout_18"));
        xmax = new QLineEdit(tab_4);
        xmax->setObjectName(QStringLiteral("xmax"));
        xmax->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        xmax->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_18->addWidget(xmax);

        ymax = new QLineEdit(tab_4);
        ymax->setObjectName(QStringLiteral("ymax"));
        ymax->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        ymax->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_18->addWidget(ymax);

        zmax = new QLineEdit(tab_4);
        zmax->setObjectName(QStringLiteral("zmax"));
        zmax->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        zmax->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_18->addWidget(zmax);


        horizontalLayout_4->addLayout(verticalLayout_18);

        verticalLayout_26 = new QVBoxLayout();
        verticalLayout_26->setSpacing(6);
        verticalLayout_26->setObjectName(QStringLiteral("verticalLayout_26"));
        label_38 = new QLabel(tab_4);
        label_38->setObjectName(QStringLiteral("label_38"));

        verticalLayout_26->addWidget(label_38);

        label_39 = new QLabel(tab_4);
        label_39->setObjectName(QStringLiteral("label_39"));

        verticalLayout_26->addWidget(label_39);

        label_43 = new QLabel(tab_4);
        label_43->setObjectName(QStringLiteral("label_43"));

        verticalLayout_26->addWidget(label_43);


        horizontalLayout_4->addLayout(verticalLayout_26);


        verticalLayout_29->addLayout(horizontalLayout_4);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        verticalLayout_21 = new QVBoxLayout();
        verticalLayout_21->setSpacing(6);
        verticalLayout_21->setObjectName(QStringLiteral("verticalLayout_21"));
        label_34 = new QLabel(tab_4);
        label_34->setObjectName(QStringLiteral("label_34"));

        verticalLayout_21->addWidget(label_34);

        label_35 = new QLabel(tab_4);
        label_35->setObjectName(QStringLiteral("label_35"));

        verticalLayout_21->addWidget(label_35);

        label_36 = new QLabel(tab_4);
        label_36->setObjectName(QStringLiteral("label_36"));

        verticalLayout_21->addWidget(label_36);


        horizontalLayout_5->addLayout(verticalLayout_21);

        verticalLayout_22 = new QVBoxLayout();
        verticalLayout_22->setSpacing(6);
        verticalLayout_22->setObjectName(QStringLiteral("verticalLayout_22"));
        xRes = new QLabel(tab_4);
        xRes->setObjectName(QStringLiteral("xRes"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(xRes->sizePolicy().hasHeightForWidth());
        xRes->setSizePolicy(sizePolicy);

        verticalLayout_22->addWidget(xRes);

        yRes = new QLabel(tab_4);
        yRes->setObjectName(QStringLiteral("yRes"));
        sizePolicy.setHeightForWidth(yRes->sizePolicy().hasHeightForWidth());
        yRes->setSizePolicy(sizePolicy);

        verticalLayout_22->addWidget(yRes);

        zRes = new QLabel(tab_4);
        zRes->setObjectName(QStringLiteral("zRes"));
        sizePolicy.setHeightForWidth(zRes->sizePolicy().hasHeightForWidth());
        zRes->setSizePolicy(sizePolicy);

        verticalLayout_22->addWidget(zRes);


        horizontalLayout_5->addLayout(verticalLayout_22);

        verticalLayout_25 = new QVBoxLayout();
        verticalLayout_25->setSpacing(6);
        verticalLayout_25->setObjectName(QStringLiteral("verticalLayout_25"));
        label_40 = new QLabel(tab_4);
        label_40->setObjectName(QStringLiteral("label_40"));

        verticalLayout_25->addWidget(label_40);

        label_41 = new QLabel(tab_4);
        label_41->setObjectName(QStringLiteral("label_41"));

        verticalLayout_25->addWidget(label_41);

        label_42 = new QLabel(tab_4);
        label_42->setObjectName(QStringLiteral("label_42"));

        verticalLayout_25->addWidget(label_42);


        horizontalLayout_5->addLayout(verticalLayout_25);


        verticalLayout_29->addLayout(horizontalLayout_5);


        horizontalLayout_9->addLayout(verticalLayout_29);


        verticalLayout_30->addLayout(horizontalLayout_9);

        verticalSpacer_4 = new QSpacerItem(20, 7, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_30->addItem(verticalSpacer_4);

        tabWidget->addTab(tab_4, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QStringLiteral("tab_2"));
        verticalLayout_27 = new QVBoxLayout(tab_2);
        verticalLayout_27->setSpacing(6);
        verticalLayout_27->setContentsMargins(11, 11, 11, 11);
        verticalLayout_27->setObjectName(QStringLiteral("verticalLayout_27"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        verticalLayout_13 = new QVBoxLayout();
        verticalLayout_13->setSpacing(6);
        verticalLayout_13->setObjectName(QStringLiteral("verticalLayout_13"));
        label = new QLabel(tab_2);
        label->setObjectName(QStringLiteral("label"));

        verticalLayout_13->addWidget(label);

        label_2 = new QLabel(tab_2);
        label_2->setObjectName(QStringLiteral("label_2"));

        verticalLayout_13->addWidget(label_2);

        label_3 = new QLabel(tab_2);
        label_3->setObjectName(QStringLiteral("label_3"));

        verticalLayout_13->addWidget(label_3);

        label_4 = new QLabel(tab_2);
        label_4->setObjectName(QStringLiteral("label_4"));

        verticalLayout_13->addWidget(label_4);


        horizontalLayout_2->addLayout(verticalLayout_13);

        verticalLayout_12 = new QVBoxLayout();
        verticalLayout_12->setSpacing(6);
        verticalLayout_12->setObjectName(QStringLiteral("verticalLayout_12"));
        NumThreads = new QLineEdit(tab_2);
        NumThreads->setObjectName(QStringLiteral("NumThreads"));
        NumThreads->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        NumThreads->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_12->addWidget(NumThreads);

        SimMode = new QLineEdit(tab_2);
        SimMode->setObjectName(QStringLiteral("SimMode"));
        SimMode->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        SimMode->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_12->addWidget(SimMode);

        SimTimeStep = new QLineEdit(tab_2);
        SimTimeStep->setObjectName(QStringLiteral("SimTimeStep"));
        SimTimeStep->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        SimTimeStep->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_12->addWidget(SimTimeStep);

        OutputFreq = new QLineEdit(tab_2);
        OutputFreq->setObjectName(QStringLiteral("OutputFreq"));
        OutputFreq->setStyleSheet(QStringLiteral("background-color: rgb(255, 255, 255);"));
        OutputFreq->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        verticalLayout_12->addWidget(OutputFreq);


        horizontalLayout_2->addLayout(verticalLayout_12);


        verticalLayout_27->addLayout(horizontalLayout_2);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_27->addItem(verticalSpacer);

        tabWidget->addTab(tab_2, QString());

        verticalLayout_10->addWidget(tabWidget);

        progressBar = new QProgressBar(centralWidget);
        progressBar->setObjectName(QStringLiteral("progressBar"));
        progressBar->setValue(0);

        verticalLayout_10->addWidget(progressBar);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_8->addItem(horizontalSpacer_2);

        RunSimulation = new QPushButton(centralWidget);
        RunSimulation->setObjectName(QStringLiteral("RunSimulation"));

        horizontalLayout_8->addWidget(RunSimulation);

        StopSim = new QPushButton(centralWidget);
        StopSim->setObjectName(QStringLiteral("StopSim"));

        horizontalLayout_8->addWidget(StopSim);


        verticalLayout_10->addLayout(horizontalLayout_8);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 369, 21));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuFile->addAction(Import);

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", nullptr));
        Import->setText(QApplication::translate("MainWindow", "Import", nullptr));
        SimNameLabel->setText(QApplication::translate("MainWindow", "Simulation Name:", nullptr));
        OutputLabel->setText(QApplication::translate("MainWindow", "Output Directory:", nullptr));
        PathLabel->setText(QApplication::translate("MainWindow", "Path File Directory: ", nullptr));
        SimName->setText(QString());
        OutputDir->setText(QString());
        PathDir->setText(QString());
        BrowseOutput->setText(QApplication::translate("MainWindow", "Browse", nullptr));
        BrowsePath->setText(QApplication::translate("MainWindow", "Browse", nullptr));
        label_37->setText(QApplication::translate("MainWindow", "Simulation Parameters", nullptr));
        label_6->setText(QApplication::translate("MainWindow", "Solidification Temperature", nullptr));
        label_7->setText(QApplication::translate("MainWindow", "Specific Heat", nullptr));
        label_8->setText(QApplication::translate("MainWindow", "Thermal Conductivity", nullptr));
        label_9->setText(QApplication::translate("MainWindow", "Density", nullptr));
        SolidTemp->setText(QApplication::translate("MainWindow", "1528.0", nullptr));
        SpecificHeat->setText(QApplication::translate("MainWindow", "600.00", nullptr));
        ThermalConductivity->setText(QApplication::translate("MainWindow", "26.6", nullptr));
        Density->setText(QApplication::translate("MainWindow", "7451.0", nullptr));
        label_11->setText(QApplication::translate("MainWindow", "K", nullptr));
        label_12->setText(QApplication::translate("MainWindow", "J/Kg-K", nullptr));
        label_13->setText(QApplication::translate("MainWindow", "W/mK", nullptr));
        label_14->setText(QApplication::translate("MainWindow", "kg/m3", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("MainWindow", "Material Prop.", nullptr));
        label_15->setText(QApplication::translate("MainWindow", "Beam Width Ax", nullptr));
        label_16->setText(QApplication::translate("MainWindow", "Beam Width Ay", nullptr));
        label_17->setText(QApplication::translate("MainWindow", "Beam Width Az", nullptr));
        label_18->setText(QApplication::translate("MainWindow", "Absorption Efficiency", nullptr));
        label_19->setText(QApplication::translate("MainWindow", "Beam Power", nullptr));
        label_5->setText(QApplication::translate("MainWindow", "Preheat Temperature", nullptr));
        BeamAx->setText(QApplication::translate("MainWindow", "10", nullptr));
        BeamAy->setText(QApplication::translate("MainWindow", "10", nullptr));
        BeamAz->setText(QApplication::translate("MainWindow", "1", nullptr));
        AbsorptionEff->setText(QApplication::translate("MainWindow", "1.0", nullptr));
        BeamPower->setText(QApplication::translate("MainWindow", "1200", nullptr));
        PreheatTemp->setText(QApplication::translate("MainWindow", "1273.0", nullptr));
        label_20->setText(QApplication::translate("MainWindow", "\302\265m", nullptr));
        label_21->setText(QApplication::translate("MainWindow", "\302\265m", nullptr));
        label_22->setText(QApplication::translate("MainWindow", "\302\265m", nullptr));
        label_23->setText(QApplication::translate("MainWindow", "%", nullptr));
        label_24->setText(QApplication::translate("MainWindow", "W", nullptr));
        label_10->setText(QApplication::translate("MainWindow", "K", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(tab_3), QApplication::translate("MainWindow", "Beam Prop.", nullptr));
        label_25->setText(QApplication::translate("MainWindow", "X min", nullptr));
        label_26->setText(QApplication::translate("MainWindow", "Y min", nullptr));
        label_27->setText(QApplication::translate("MainWindow", "Z min", nullptr));
        label_31->setText(QApplication::translate("MainWindow", "imax", nullptr));
        label_32->setText(QApplication::translate("MainWindow", "jmax", nullptr));
        label_33->setText(QApplication::translate("MainWindow", "kmax", nullptr));
        xmin->setText(QApplication::translate("MainWindow", "-1", nullptr));
        ymin->setText(QApplication::translate("MainWindow", "-1", nullptr));
        zmin->setText(QApplication::translate("MainWindow", "-.5", nullptr));
        imax->setText(QApplication::translate("MainWindow", "280", nullptr));
        jmax->setText(QApplication::translate("MainWindow", "150", nullptr));
        kmax->setText(QApplication::translate("MainWindow", "1", nullptr));
        label_44->setText(QApplication::translate("MainWindow", "mm", nullptr));
        label_45->setText(QApplication::translate("MainWindow", "mm", nullptr));
        label_46->setText(QApplication::translate("MainWindow", "mm", nullptr));
        label_47->setText(QApplication::translate("MainWindow", "Pts", nullptr));
        label_48->setText(QApplication::translate("MainWindow", "Pts", nullptr));
        label_49->setText(QApplication::translate("MainWindow", "Pts", nullptr));
        label_28->setText(QApplication::translate("MainWindow", "X max", nullptr));
        label_29->setText(QApplication::translate("MainWindow", "Y max", nullptr));
        label_30->setText(QApplication::translate("MainWindow", "Z max", nullptr));
        xmax->setText(QApplication::translate("MainWindow", "6", nullptr));
        ymax->setText(QApplication::translate("MainWindow", "2", nullptr));
        zmax->setText(QApplication::translate("MainWindow", "0.0", nullptr));
        label_38->setText(QApplication::translate("MainWindow", "mm", nullptr));
        label_39->setText(QApplication::translate("MainWindow", "mm", nullptr));
        label_43->setText(QApplication::translate("MainWindow", "mm", nullptr));
        label_34->setText(QApplication::translate("MainWindow", "X Resolution", nullptr));
        label_35->setText(QApplication::translate("MainWindow", "Y Resolution", nullptr));
        label_36->setText(QApplication::translate("MainWindow", "Z Resolution", nullptr));
        xRes->setText(QString());
        yRes->setText(QString());
        zRes->setText(QString());
        label_40->setText(QApplication::translate("MainWindow", "\302\265m", nullptr));
        label_41->setText(QApplication::translate("MainWindow", "\302\265m", nullptr));
        label_42->setText(QApplication::translate("MainWindow", "\302\265m", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(tab_4), QApplication::translate("MainWindow", "Simulation Bounds", nullptr));
        label->setText(QApplication::translate("MainWindow", "Number of Threads", nullptr));
        label_2->setText(QApplication::translate("MainWindow", "Simulation Mode", nullptr));
        label_3->setText(QApplication::translate("MainWindow", "Simulation Timestep", nullptr));
        label_4->setText(QApplication::translate("MainWindow", "File Output Frequency", nullptr));
        NumThreads->setText(QApplication::translate("MainWindow", "4", nullptr));
        SimMode->setText(QApplication::translate("MainWindow", "2", nullptr));
        SimTimeStep->setText(QApplication::translate("MainWindow", ".0001", nullptr));
        OutputFreq->setText(QApplication::translate("MainWindow", "10", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QApplication::translate("MainWindow", "Settings", nullptr));
        RunSimulation->setText(QApplication::translate("MainWindow", "Start Sim", nullptr));
        StopSim->setText(QApplication::translate("MainWindow", "Stop Sim", nullptr));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
