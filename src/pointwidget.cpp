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

#include "pointwidget.h"
#include "ui_pointwidget.h"

#include "Data.h"
#include "truncatedialog.h"
#include "bsplineinterpolation.h"
#include "settings.h"
#include "measurementpoint.h"

#include <QAbstractItemModel>
#include <QList>
#include <QFileInfo>
#include <QEvent>
#include <QFileDialog>
#include <QTableWidgetItem>
#include <QMessageBox>

#include <cmath>
#include <vector>
#include <limits>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


#undef min
#undef max


class PointModel : public QAbstractTableModel {
public:
	typedef QList<Wm3::Vector3<double> > PointDataType;

	PointModel() {}

	int rowCount(const QModelIndex& parent = QModelIndex()) const {
		return m_data.size();
	}

	int columnCount(const QModelIndex& parent = QModelIndex()) const {
		return 3;
	}

	QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const {

        if( role == Qt::TextAlignmentRole )
            return Qt::AlignCenter;

		if (!index.isValid() || (role != Qt::DisplayRole) || index.row() >= m_data.size() || index.row() < 0 || index.column() < 0 || index.column() > 2)
			return QVariant();

		// Return the data to which index points.
		return m_data[index.row()][index.column()];
	}

	QVariant headerData(int section, Qt::Orientation orientation, int role) const {
		if (role != Qt::DisplayRole)
			return QVariant();

		if (orientation == Qt::Horizontal) {
			switch (section) {
				case 0:
					return tr("X");
				case 1:
					return tr("Y");
				case 2:
					return tr("Z");
				default:
					return QVariant();
			}
		}

		return QVariant();
	}

	PointDataType getAllData() const {
		return m_data;
	}

	bool setData(const QModelIndex& index, const QVariant& value, int role=Qt::EditRole) {
		if (index.isValid() && (role == Qt::EditRole) && (index.row() >= 0) && (index.row() < m_data.size()) && (index.column() >= 0) && (index.column() < 3)) {
			int row = index.row();

			Wm3::Vector3<double> v = m_data.value(row);
			v[index.column()] = value.toDouble();

			m_data.replace(row, v);
			emit(dataChanged(index, index));

			return true;
		}

		return false;
	}

	bool insertRows(int position, int rows, const QModelIndex& index = QModelIndex()) {
		beginInsertRows(QModelIndex(), position, position+rows-1);

		for (int row=0; row < rows; row++)
			m_data.insert(position, Wm3::Vector3<double>::ZERO);

		endInsertRows();
		
		return true;
	}

	bool removeRows(int position, int rows, const QModelIndex& index = QModelIndex()) {
		if (rows == 0)
			return true;

		if (m_data.size() < rows)
			return false;

		beginRemoveRows(QModelIndex(), position, position+rows-1);

		for (int row=0; row < rows; ++row)
			m_data.removeAt(position);

		endRemoveRows();

		return true;
	}

	void clear() {
        removeRows(0, m_data.size());
	}

private:
	 PointDataType m_data;
};



PointWidget::PointWidget(QWidget* parent) :
		QDialog(parent),
    m_settings(Settings::getSettings()),
		ui(new Ui::PointWidget)
{
    ui->setupUi(this);
    connect(ui->b_clear,   SIGNAL(clicked()), this, SLOT(clearTable()));
    connect(ui->b_load,    SIGNAL(clicked()), this, SLOT(loadControlPoints()) );
    connect(ui->b_save,    SIGNAL(clicked()), this, SLOT(saveControlPoints()) );
	connect(ui->b_truncate, SIGNAL(clicked()), this, SLOT(truncateSpline()));
    
	PointModel* model = new PointModel();
	ui->pointTable->setModel(model);
	
    m_data    = 0;
    m_bSpline = 0;

	clearTable();

    setWindowTitle("Control Points");

	m_lastDir = m_settings.getDefaultDir();
}

PointWidget::~PointWidget()
{
	PointModel* model = dynamic_cast<PointModel*>(ui->pointTable->model());
    delete ui;
	delete model;
}

void PointWidget::reset()
{
	clearTable();

	m_bSpline = NULL;
	m_data = NULL;

	ui->b_load->setEnabled(false);
}

void PointWidget::setInterpolator(BSplineInterpolation* bSpline)
{
	m_bSpline = bSpline;

	if (!m_bSpline)
		return;

	ui->b_load->setEnabled(true);

	PointModel* model = dynamic_cast<PointModel*>(ui->pointTable->model());
	if (model->rowCount() == 0)
		return;

	ui->b_clear->setEnabled(true);
	ui->b_save->setEnabled(true);
	ui->b_truncate->setEnabled(true);
}

std::vector< Wm3::Vector3<double> > PointWidget::getControlPoints()
{
    // TODO This is ugly!!!!!!!!
	PointModel* model = dynamic_cast<PointModel*>(ui->pointTable->model());
	return model->getAllData().toVector().toStdVector();
}

void PointWidget::clearTable()
{
	PointModel* model = dynamic_cast<PointModel*>(ui->pointTable->model());
    model->clear();

    if (m_bSpline)
		m_bSpline->reset();

	ui->b_clear->setEnabled(false);
	ui->b_save->setEnabled(false);
	ui->b_truncate->setEnabled(false);

    emit pointsChanged();
}

void PointWidget::addControlPoint(const Wm3::Vector3<double>& p)
{
	if (!m_bSpline || !m_data)
		return;

	PointModel* model = dynamic_cast<PointModel*>(ui->pointTable->model());
	
	int row = model->rowCount();
	model->insertRow(row, QModelIndex());
	for (int i = 0; i < 3; ++i) {
		QModelIndex index = model->index(row, i, QModelIndex());
		model->setData(index, p[i], Qt::EditRole);
	}

	ui->b_clear->setEnabled(true);
	ui->b_save->setEnabled(true);
	ui->b_truncate->setEnabled(true);
}

void PointWidget::loadControlPoints()
{
    if (!m_bSpline || !m_data)
		return;

	QString fileName = QFileDialog::getOpenFileName(this,
             tr("Load Control Points"), m_lastDir,
             tr("Segmentation Files (*.txt *.xml)"));

	if (fileName.isEmpty())
		return;

	clearTable();

	if (m_settings.getRememberLastDir()) {
		QFileInfo info(fileName);
		m_lastDir = info.path();
	}

	// Try to determine file type
	int dotPos = fileName.lastIndexOf(".");
	QString ext;
	if (dotPos == -1)
		ext = ".txt";
	else {
		ext = fileName.right(fileName.length() - dotPos - 1);
	}
	
	if (ext.compare("xml", Qt::CaseInsensitive) != 0) {
		// Try to load point text file
		
		// Parse spline file
		QByteArray fn = fileName.toLocal8Bit();
		std::ifstream file(fn.constData());
		if (!file.is_open() || !file.good()) {
			QMessageBox::critical(this, "Error", "Could not open spline text file");
			return;
		}

		Wm3::Vector3<double> pt;
		std::stringstream line;
		std::string tmp;

		do {
			if (!std::getline(file, tmp))
				break;

			line.clear();
			line << tmp.c_str();

			line >> pt[0] >> pt[1] >> pt[2];
			if (line.bad())
				break;

			addControlPoint(pt);

		} while (!file.eof() && file.good());

		file.close();
	} else {
		QByteArray fn = fileName.toLocal8Bit();
		std::ifstream file(fn.constData());
		if (!file.is_open() || !file.good()) {
			QMessageBox::critical(this, "Error", "Could not open spline XML file");
			return;
		}

		while (!file.eof() && file.good()) {
			std::string line;
			std::getline(file, line);
			if (line.empty())
				continue;

			size_t begin = line.find("<pos>");
			if (begin == std::string::npos)
				continue;

			size_t end = line.find("</pos>");
			if (end == std::string::npos) {
				std::cerr << "Encountered invalid line " << line << std::endl;
				continue;
			}

			begin += 5;
			end -= 7;

			std::istringstream parser(line.substr(begin, end - begin + 1));

			Wm3::Vector3<double> pt;
			parser >> pt[0] >> pt[1] >> pt[2];

			if (parser.bad()) {
				std::cerr << "Error parsing string " << line << std::endl;
				continue;
			}

			addControlPoint(pt);
		}

		file.close();
	}

	ui->b_clear->setEnabled(true);
	ui->b_save->setEnabled(true);
	ui->b_truncate->setEnabled(true);

	emit pointsChanged();
}

void PointWidget::saveControlPoints()
{
    QString fileName = QFileDialog::getSaveFileName(this,
			 tr("Save Control Points"), m_lastDir,
             tr("Text files (*.txt)"));

	if (fileName.isEmpty())
		return;

	if (m_settings.getRememberLastDir()) {
		QFileInfo info(fileName);
		m_lastDir = info.path();
	}

    // Parse spline file
	QByteArray fn = fileName.toLocal8Bit();
	std::ofstream file(fn.constData());
    if (!file.good()) {
        std::cerr << "Could not open spline file" << std::endl;
        return;
    }

	PointModel* model = dynamic_cast<PointModel*>(ui->pointTable->model());

	file.precision(6);
	for (PointModel::PointDataType::const_iterator i = model->getAllData().constBegin();
		i != model->getAllData().constEnd();
		++i) {
			file << std::fixed
		     << (*i)[0] << " "
             << (*i)[1] << " "
             << (*i)[2] << std::endl;
    }
    file.close();
}

void PointWidget::improveControlPoints()
{
    if( !m_bSpline || !m_data)
        return;

    if( m_bSpline->getLength() < 0 )
        return;

	double length = m_bSpline->getLength() * m_data->getVoxelSize();
    double evalInterval = m_settings.getSharpnessEvalInterval();
    double delta = evalInterval / length;
    std::vector<MeasurementPoint*> measurements;
    for(double i=0.0; i<=1.0; i += delta)
    {
        MeasurementPoint* mp = new MeasurementPoint(i, m_bSpline, m_data);
        mp->improvePoint();
        measurements.push_back(mp);
    }

    Wm3::Vector3<double> p;
    for(size_t i=0; i<measurements.size(); i++)
    {
        MeasurementPoint* mp = measurements.at(i);
        if( i > 0 && i < measurements.size()-1 )
        {
            double dist = (p - mp->getPositionInVolume()).Length() * m_data->getVoxelSize();
            if(dist > 1.5 * evalInterval)
            {
                MeasurementPoint* prevP = measurements.at(i-1);
                MeasurementPoint* nextP = measurements.at(i+1);
                mp->setPositionInVolume( 0.5*(prevP->getPositionInVolume()+nextP->getPositionInVolume()) );
            }
        }
        p = mp->getPositionInVolume();
    }

    clearTable();
     for(size_t i=0; i<measurements.size(); i++)
         addControlPoint(measurements.at(i)->getPositionInVolume());

    for(size_t i=0; i<measurements.size(); i++)
        delete measurements.at(i);
    measurements.clear();
    emit pointsChanged();

    QMessageBox::information(NULL, "Points improved", "Control Points have been successfully improved!");
}

void PointWidget::resizeEvent(QResizeEvent *event)
{
    QWidget::resizeEvent(event);
	int width = ceil(double( ui->pointTable->width() - 20 ) / 3.0);
    for(int i=0; i < 3; i++)
        ui->pointTable->setColumnWidth(i, width);
}

void PointWidget::truncateSpline()
{
	if (!m_bSpline || !m_data)
		return;

	double length = m_bSpline->getLength();

	TruncateDialog dialog(this);

	dialog.setRange(0, length * m_data->getVoxelSize());

	if (dialog.exec() == QDialog::Rejected)
		return;

	double start = dialog.getStart() / m_data->getVoxelSize();
	double end = dialog.getEnd() / m_data->getVoxelSize();

	// No points left
	if (start == end) {
		clearTable();
		return;
	}

	// Full length
	if ((start < dialog.getStep() / m_data->getVoxelSize()) && (std::abs(end - length) < dialog.getStep() / m_data->getVoxelSize())) {
		return;
	}

	// Normalise to spline index
	start /= length;
	end /= length;

	// Find closest interpolated points
	Wm3::Vector3<double> pStart = m_bSpline->pointAt(start);
	Wm3::Vector3<double> pEnd = m_bSpline->pointAt(end);

	std::vector<Wm3::Vector3<double> > points = m_bSpline->getInterpolatedPoints();
	std::vector<Wm3::Vector3<double> >::const_iterator it, itStart, itEnd;
	double distStart = std::numeric_limits<double>::max();
	double distEnd = std::numeric_limits<double>::max();
	
	for (it = points.begin(); it != points.end(); ++it) {
		double d = (pStart - *it).Length();
		if (d < distStart) {
			itStart = it;
			distStart = d;
		}

		d = (pEnd - *it).Length();
		if (d < distEnd) {
			itEnd = it;
			distEnd = d;
		}
	}

	// Change spline points to truncated points
	clearTable();
	for (it = itStart; it != itEnd + 1; ++it)
		addControlPoint(*it);

	emit pointsChanged();
}

