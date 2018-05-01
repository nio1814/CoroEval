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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>


// Forward declarations
class Settings;
class BSplineInterpolation;
class Data;
class EvaluationDialog;
class NormalPlot;
class MeasurementPoint;
class QKeyEvent;

namespace Ui {
class MainWindow;
}


class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget* parent = 0);
    ~MainWindow();

protected slots:
    void loadFile();
	void clearData();
    void setData(Data* data);
    void updateMeasurementPoints();
    void updatePosition(size_t p) const;
	void recalculate(bool recache = true);
	void measurementsEdit(size_t p);

	void showNormalPlot();
	void showEvaluationDialog();

	void refineSeg();
	void optCurPoint();
	void optAllPoints();
    void exportSegmentationTriangulation();
    void exportTubeTriangulation();
    void exportCenterline();
    void exportSurfacePoints();

	void configure();
	void loadConfig();
	void saveConfig();

	void about();
	void aboutQt();

protected:
	virtual void keyPressEvent(QKeyEvent* event);

private:
	Settings&						m_settings;
    Ui::MainWindow*					ui;

    BSplineInterpolation*           m_bSpline;
    Data*                           m_data;

    EvaluationDialog*               m_evalDialog;
    NormalPlot*                     m_normalPlot;
    std::vector<MeasurementPoint*>  m_measurements;

	QString							m_lastLoadDir;

    void reset(bool resetSpline = true);
};

#endif // MAINWINDOW_H

