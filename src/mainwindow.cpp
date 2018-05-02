/*
 * CoroEval
 *
 * An evaluation tool for coronary artery reconstructions.
 *
 * Copyright © 2014:
 *
 * Christoph Forman, Universität Erlangen-Nürnberg
 * christoph.forman@cs.fau.de
 *
 * Chris Schwemmer, Universität Erlangen-Nürnberg
 * chris.schwemmer@cs.fau.de
 *
 * Jens Wetzl, Universität Erlangen-Nürnberg
 * jens.wetzl@cs.fau.de
 *
 * CoroEval is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CoroEval is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CoroEval.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "settings.h"
#include "settingsdialog.h"
#include "bsplineinterpolation.h"
#include "Data.h"
#include "evaluationdialog.h"
#include "normalplot.h"
#include "measurementpoint.h"
#include "loaddialog.h"
#include "os.h"
#include "MPFilter.h"

#ifdef ENABLE_MESH_EXPORT
#include "Triangulation.h"
#endif

#include <CImg.h>

#include <QTime>
#include <QFile>
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QKeyEvent>

#include <iostream>
#include <fstream>
#include <algorithm>


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    m_settings(Settings::getSettings()),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Get StyleSheet from file :)
    QFile file(":/res/style.qss");
    file.open(QFile::ReadOnly);
    QString styleSheet = QLatin1String(file.readAll());
    setStyleSheet(styleSheet);

    ui->viewer->setStyleSheet(styleSheet);

	// Try to load default config if available
	m_settings.load("config.txt");

    m_data = new Data();
    m_bSpline = new BSplineInterpolation();

    m_evalDialog = new EvaluationDialog(this);
    m_evalDialog->setInterpolator(m_bSpline);

	if (m_settings.getShowNormals()) {
		m_normalPlot = new NormalPlot(this);
		m_normalPlot->setStyleSheet(styleSheet);
		m_normalPlot->show();
	} else
		m_normalPlot = NULL;

    ui->plot2->setMode(SharpnessPlot::Diameter);

    connect(ui->profile, SIGNAL(showSlice(const Wm3::Vector3<double>&)), ui->viewer, SLOT(showSlice(const Wm3::Vector3<double>&)));
    connect(ui->viewer, SIGNAL(updateProfile()), this, SLOT(updateMeasurementPoints()));
    connect(ui->profile, SIGNAL(positionChanged(size_t)), this, SLOT(updatePosition(size_t)));
	connect(ui->profile, SIGNAL(measurementsEdit(size_t)), this, SLOT(measurementsEdit(size_t)));

	connect(ui->actionOpen, SIGNAL(triggered()), this, SLOT(loadFile()));
	connect(ui->actionExit, SIGNAL(triggered()), this, SLOT(close()));
	connect(ui->actionPoints, SIGNAL(triggered()), ui->viewer, SLOT(showPointsDialog()));
	connect(ui->actionNormals, SIGNAL(triggered()), this, SLOT(showNormalPlot()));
	connect(ui->actionEvaluationDialog, SIGNAL(triggered()), this, SLOT(showEvaluationDialog()));
	connect(ui->actionRefine_Segmentation_in_Plane, SIGNAL(triggered()), this, SLOT(refineSeg()));
	connect(ui->actionOptimise_Current_Point, SIGNAL(triggered()), this, SLOT(optCurPoint()));
	connect(ui->actionOptimise_All_Points, SIGNAL(triggered()), this, SLOT(optAllPoints()));
    connect(ui->actionExport_Segmentation_Mesh, SIGNAL(triggered(bool)), this, SLOT(exportSegmentationTriangulation()));
    connect(ui->actionExport_Tube_Mesh, SIGNAL(triggered(bool)), this, SLOT(exportTubeTriangulation()));
    connect(ui->actionExport_Centerline, SIGNAL(triggered(bool)), this, SLOT(exportCenterline()));
    connect(ui->actionExport_Surface_Points, SIGNAL(triggered(bool)), this, SLOT(exportSurfacePoints()));
	connect(ui->actionConfigure, SIGNAL(triggered()), this, SLOT(configure()));
	connect(ui->actionLoad, SIGNAL(triggered()), this, SLOT(loadConfig()));
	connect(ui->actionSave, SIGNAL(triggered()), this, SLOT(saveConfig()));
	connect(ui->actionCoroEval, SIGNAL(triggered()), this, SLOT(about()));
	connect(ui->actionQt, SIGNAL(triggered()), this, SLOT(aboutQt()));

#ifndef ENABLE_MESH_EXPORT
    ui->menuExport_Triangulation->setVisible(false);
#endif

	m_lastLoadDir = m_settings.getDefaultDir();
}

MainWindow::~MainWindow()
{
	if (m_normalPlot)
		delete m_normalPlot;
	delete m_evalDialog;
	delete m_bSpline;
	delete m_data;

    delete ui;
}

void MainWindow::loadFile()
{
    LoadDialog dialog(this);
	connect(&dialog, SIGNAL(clearData()), this, SLOT(clearData()) );
    connect(&dialog, SIGNAL(sendData(Data*)), this, SLOT(setData(Data*)) );
    connect(&dialog, SIGNAL(sendFileName(QString)), this, SLOT(setWindowTitle(QString)) );
	dialog.setLastLoadDir(m_lastLoadDir);

    dialog.exec();

	disconnect(this, SLOT(clearData()));
    disconnect(this, SLOT(setData(Data*)));
    disconnect(this, SLOT(setWindowTitle(QString)));
	m_lastLoadDir = dialog.getLastLoadDir();
}

void MainWindow::clearData()
{
	if (m_data->data0()) {
		ui->viewer->reset();
		ui->profile->reset();
		m_evalDialog->reset();

		delete m_data->data0();
		m_data->setData(NULL);
    }

	reset();
}

void MainWindow::setData(Data* data)
{
    // Old dataset?
    if (m_data->data0()) {
		ui->viewer->reset();
		ui->profile->reset();
		m_evalDialog->reset();

		delete m_data->data0();
    }

	m_data->clone(data);
    m_data->setVoxelSize(data->getVoxelSize());

	delete data;
	data = NULL;

    ui->viewer->setInterpolator(m_bSpline);
    ui->viewer->setData(m_data);

    ui->profile->setData(m_data);
    ui->profile->setProfileSize( 15 );

	m_evalDialog->setInterpolator(m_bSpline);
    m_evalDialog->setData(m_data);

    this->resize(2.5*ui->viewer->width(), this->height());

    reset();
}

void MainWindow::updateMeasurementPoints()
{
    if( !m_bSpline || !m_data)
        return;

    reset(false);

    if( m_bSpline->getLength() < 0 )
        return;

    QTime t;
    t.start();

    double length = m_bSpline->getLength() * m_data->getVoxelSize();
    double evalInterval = m_settings.getSharpnessEvalInterval();
	size_t numSamplePoints = static_cast<size_t>(length / evalInterval);

	// Create new measurements
    for(size_t i = 0; i < numSamplePoints; ++i)
    {
        MeasurementPoint* mp = new MeasurementPoint(double(i) / double(numSamplePoints - 1), m_bSpline, m_data);
        m_measurements.push_back( mp );
    }

	// Remove 0s (Problem does not seem to appear anymore?)
	//MPFilter::purge0(m_measurements);

	// Median-filter diameters
	MPFilter::median(m_measurements, 3);
	//MPFilter::gauss(m_measurements, 3, 0.5);

	// Create lengthwise profile of the coronary
	size_t width = (static_cast<double>(m_settings.getSizeNormal()) / m_data->getVoxelSize()) * 2 + 1;

	cimg_library::CImg<float> cImage(width, m_measurements.size());
	cImage.fill(0);
	float* image = cImage.data();

	for (size_t i = 0; i < m_measurements.size(); i++) {
		const MeasurementPoint* mp = m_measurements[i];

		size_t thr = static_cast<size_t>(std::floor(mp->getDiameter() * 2 / m_data->getVoxelSize()));
		size_t pos = i * width;

		for (size_t j = 0; j < m_settings.getNumNormalsAtPoint(); ++j) {
			const std::vector<float>& n = mp->getNormal(j);

			size_t c = n.size() / 2;
			assert(c > thr);

			for (size_t k = c - thr; k < c + thr; ++k) {
				image[pos + k] += n[k] / n.size();
			}
		}
    }

	cImage.rotate(-90, 0, 0);

	// Send measurements and update everything
    ui->plot->setMeasurements(m_measurements);

    ui->plot2->setMeasurements(m_measurements);

	if (m_normalPlot)
		m_normalPlot->setMeasurements(m_measurements);

    ui->profile->setMeasurements(m_measurements, cImage, true);

	m_evalDialog->setMeasurements(m_measurements);
	m_evalDialog->updateMax();

	ui->actionRefine_Segmentation_in_Plane->setEnabled(true);
	ui->actionOptimise_Current_Point->setEnabled(true);
	ui->actionOptimise_All_Points->setEnabled(true);
    ui->menuExport_Triangulation->setEnabled(true);

    qDebug() << "done: " << QTime().addMSecs(t.elapsed()).toString("mm.ss,zzz");
}

void MainWindow::updatePosition(size_t p) const
{
    //qDebug() << "called updatePosition: " << p;
    ui->plot->setPosition(p);
    ui->plot2->setPosition(p);
	if (m_normalPlot)
		m_normalPlot->setPosition(p);
}

void MainWindow::reset(bool resetSpline)
{
    // Clear existing measurements
    for(size_t i=0; i<m_measurements.size(); i++)
        delete m_measurements.at(i);
    m_measurements.clear();

    ui->profile->setMeasurements(m_measurements);

    ui->plot->setMeasurements(m_measurements);
    ui->plot->setPosition(0);

    ui->plot2->setMeasurements(m_measurements);
    ui->plot2->setPosition(0);

	if (m_normalPlot)
		m_normalPlot->setMeasurements(m_measurements);

	m_evalDialog->setMeasurements(m_measurements);
	m_evalDialog->updateMax();

    if( m_bSpline && resetSpline )
        m_bSpline->reset();

	ui->actionRefine_Segmentation_in_Plane->setEnabled(false);
	ui->actionOptimise_Current_Point->setEnabled(false);
	ui->actionOptimise_All_Points->setEnabled(false);
    ui->menuExport_Triangulation->setEnabled(false);
}

void MainWindow::recalculate(bool recache) {
	size_t width = (static_cast<double>(m_settings.getSizeNormal()) / m_data->getVoxelSize()) * 2 + 1;

	cimg_library::CImg<float> cImage(width, m_measurements.size());
	cImage.fill(0);
	float* image = cImage.data();

	// Re-calculate mesurements
	for (size_t i = 0; i < m_measurements.size(); i++)
	{
        m_measurements[i]->recalculate();
	}

	// Remove 0s (Problem does not seem to appear anymore?)
	//MPFilter::purge0(m_measurements);

	// Median-filter measurements
	MPFilter::median(m_measurements, 3);

	// Create lengthwise profile of the coronary
	for (size_t i = 0; i < m_measurements.size(); i++) {
		MeasurementPoint* mp = m_measurements[i];

		size_t thr = static_cast<size_t>(std::floor(mp->getDiameter() * 2 / m_data->getVoxelSize()));
		size_t pos = i * width;

		for (size_t j = 0; j < m_settings.getNumNormalsAtPoint(); ++j) {
			const std::vector<float>& n = mp->getNormal(j);

			size_t c = n.size() / 2;
			assert(c > thr);

			for (size_t k = c - thr; k < c + thr; ++k) {
				image[pos + k] += n[k] / n.size();
			}
		}
    }

	cImage.rotate(-90, 0, 0);

	// Send measurements and update everything
    ui->plot->setMeasurements(m_measurements);

    ui->plot2->setMeasurements(m_measurements);

	if (m_normalPlot)
		m_normalPlot->setMeasurements(m_measurements);

    ui->profile->setMeasurements(m_measurements, cImage, recache);

	m_evalDialog->setMeasurements(m_measurements);
	m_evalDialog->updateMax();
}

void MainWindow::measurementsEdit(size_t p) {
	// For now, just recalculate everything
	recalculate(false);

	ui->profile->updateProfile(p);
}

void MainWindow::showNormalPlot() {
	if (!m_normalPlot) {
		m_normalPlot = new NormalPlot(this);
		m_normalPlot->setStyleSheet(styleSheet());
		m_normalPlot->setMeasurements(m_measurements);
		m_normalPlot->setPosition(ui->profile->getPosition());
	}

	m_normalPlot->show();
	m_normalPlot->activateWindow();
	m_normalPlot->raise();
}

void MainWindow::showEvaluationDialog() {
	m_evalDialog->show();
	m_evalDialog->activateWindow();
	m_evalDialog->raise();
}

void MainWindow::refineSeg() {
	if (!m_measurements.empty())
		ui->profile->setRefineMode(!ui->profile->getRefineMode());
}

void MainWindow::optCurPoint() {
	if (!m_measurements.empty()) {
		m_measurements[ui->profile->getPosition()]->improvePoint2d(0.1, 50);
		size_t pos = ui->profile->getPosition();
		recalculate(true);
		ui->profile->updateProfile(pos);
	}
}

void MainWindow::optAllPoints() {
	if (!m_measurements.empty()) {
		size_t pos = ui->profile->getPosition();
		for (size_t i = 0; i < m_measurements.size(); ++i)
			m_measurements[i]->improvePoint2d(0.1, 50);
		recalculate(true);
		ui->profile->updateProfile(pos);
	}
}

void MainWindow::exportSegmentationTriangulation()
{
#ifdef ENABLE_MESH_EXPORT
    if (m_measurements.empty())
        return;

    Triangulation tri(m_measurements);

    QString fileName = QFileDialog::getSaveFileName(this,
        tr("Save Mesh"), m_settings.getDefaultDir(),
        tr("OFF format (*.off);;OBJ format (*.obj);;STL format (*.stl)"));
    if (fileName.compare("") != 0)
        if (!tri.saveMesh(fileName.toStdString()))
            QMessageBox::information(this, "Export failed", QString("Mesh could not be saved to '%1'.").arg(fileName));
#endif
}

void MainWindow::exportTubeTriangulation()
{
#ifdef ENABLE_MESH_EXPORT
    if (m_measurements.empty())
        return;

    Triangulation tri(m_measurements);

    QString fileName = QFileDialog::getSaveFileName(this,
        tr("Save Mesh"), m_settings.getDefaultDir(),
        tr("OFF format (*.off);;OBJ format (*.obj);;STL format (*.stl)"));
    if (fileName.compare("") != 0)
        if (!tri.saveTubeMesh(fileName.toStdString()))
            QMessageBox::information(this, "Export failed", QString("Mesh could not be saved to '%1'.").arg(fileName));
#endif
}

void MainWindow::exportCenterline()
{
    if (m_measurements.empty())
        return;

    QString fileName = QFileDialog::getSaveFileName(this,
        tr("Save Centerline"), m_settings.getDefaultDir(),
        tr("Text (*.txt)"));
    if (fileName.compare("") != 0)
    {
        std::ofstream outputStream(fileName.toStdString().c_str());

        if (!outputStream.good())
        {
            QMessageBox::information(this, "Export failed", QString("Centerline could not be saved to '%1'.").arg(fileName));
            outputStream.close();
            return;
        }

        for (size_t i = 0; i < m_measurements.size(); ++i)
        {
            outputStream << m_measurements[i]->getImprovedPositionInVolume() << std::endl;
        }

        outputStream.close();
    }
}

void MainWindow::exportSurfacePoints()
{
#ifdef ENABLE_MESH_EXPORT
    if (m_measurements.empty())
        return;

    Triangulation tri(m_measurements);

    QString fileName = QFileDialog::getSaveFileName(this,
        tr("Save Surface Points"), m_settings.getDefaultDir(),
        tr("Text (*.txt)"));
    if (fileName.compare("") != 0)
        if (!tri.saveSurfacePoints(fileName.toStdString()))
            QMessageBox::information(this, "Export failed", QString("Surface points could not be saved to '%1'.").arg(fileName));
#endif
}

void MainWindow::configure() {
	SettingsDialog dialog(this);

	dialog.setupFromConfig();

	if (dialog.exec() == QDialog::Rejected)
		return;

	double e = m_settings.getSharpnessEvalInterval();

	dialog.updateConfig();

	if (fabs(e - m_settings.getSharpnessEvalInterval()) < 1e-7)
		recalculate();
	else
		updateMeasurementPoints();
}

void MainWindow::loadConfig() {
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Load Config File"), m_settings.getDefaultDir(),
		tr("Config Files (*.txt)"));

	if (fileName.isEmpty())
		return;

	if (!m_settings.load(fileName))
		QMessageBox::critical(this, "Error", "Could not load config file");
}

void MainWindow::saveConfig() {
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save Config File"), m_settings.getDefaultDir(),
		tr("Config Files (*.txt)"));

	if (fileName.isEmpty())
		return;

	if (!m_settings.save(fileName))
		QMessageBox::critical(this, "Error", "Could not save config file");
}

void MainWindow::keyPressEvent(QKeyEvent* event) {
	QMainWindow::keyPressEvent(event);
}

void MainWindow::about() {
	QString msg;
	msg = "CoroEval\n\nA multi-platform, multi-modality evaluation tool for coronary reconstructions.\n\n";
	msg += "Copyright ©2014, Christoph Forman, Chris Schwemmer & Jens Wetzl.\nPattern Recognition Lab, Universität Erlangen-Nürnberg, Germany.\n\n";
	msg += "This software is licensed under the GNU General Public License. For more details see license.txt.\n\n";
    msg += "If you use CoroEval in your research, please cite its companion publication as follows:\nCoroEval: a multi-platform, multi-modality tool for the evaluation of 3D coronary vessel reconstructions.\nC Schwemmer, C Forman, J Wetzl, A Maier and J Hornegger.\nPhys. Med. Biol. 59 (2014) 5163-5174.";

	QMessageBox::about(this, "About CoroEval", msg);
}

void MainWindow::aboutQt() {
	QMessageBox::aboutQt(this, "About Qt");
}
