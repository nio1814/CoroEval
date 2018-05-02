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

#ifndef POINTWIDGET_H
#define POINTWIDGET_H

#include "Wm3Vector3.h"

#include <QDialog>
#include <QWidget>
#include <QString>

#include <vector>

// Forward declarations
class Data;
class BSplineInterpolation;
class QResizeEvent;
class Settings;

namespace Ui {
class PointWidget;
}


class PointWidget : public QDialog
{
    Q_OBJECT
    
public:
    explicit PointWidget(QWidget* parent = 0);
    ~PointWidget();

	void reset();

    void addControlPoint(const Wm3::Vector3<double>& p);
    std::vector<Wm3::Vector3<double> > getControlPoints();

    void setData(Data* data) { m_data = data; clearTable(); }
    void setInterpolator(BSplineInterpolation* bSpline);

public slots:
    void clearTable();
    void loadControlPoints();
    void saveControlPoints();
    void improveControlPoints();
	void truncateSpline();

signals:
    void pointsChanged();

protected:
    void resizeEvent(QResizeEvent* event);

private:
	Settings& m_settings;
    Ui::PointWidget* ui;
    Data* m_data;
    BSplineInterpolation* m_bSpline;
	QString m_lastDir;
};

#endif // POINTWIDGET_H

