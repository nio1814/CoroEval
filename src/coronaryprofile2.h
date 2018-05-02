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

#ifndef CORONARYPROFILE2_H
#define CORONARYPROFILE2_H

#include <algorithm>

#include <QWidget>

#include "Wm3Vector3.h"
#include "Wm3Quaternion.h"

#define cimg_display 0
#include "CImg.h"

#include <vector>
#include <utility>


// Forward declarations
class Data;
class MeasurementPoint;
class Settings;

namespace Ui {
class CoronaryProfile2;
}

class CoronaryProfile2 : public QWidget
{
    Q_OBJECT
    
public:
    explicit CoronaryProfile2(QWidget* parent = 0);
    ~CoronaryProfile2();

	void reset();

    void setData(Data* data);
    void setProfileSize(size_t size);

	void setMeasurements(const std::vector<MeasurementPoint*>& m);

	size_t getPosition() const { return m_pos; }

public slots:
    void updateProfile(int pos);
    void update(bool recache);

signals:
    void positionChanged(size_t) const;
    void showSlice(const Wm3::Vector3<double>&) const;

private:
	const Settings& m_settings;
    Ui::CoronaryProfile2* ui;

    QImage* m_overlayImage;
    Data*   m_profile;

    size_t m_pos;
    size_t m_length;
    Data* m_data;

	std::vector<cimg_library::CImg<float> > m_planeCache;
	std::vector<std::pair<float, float> > m_scaleCache;

    std::vector<MeasurementPoint*> m_measurements;
	
	void drawAllPlanes();
    void drawPlane();
    void drawOverlay();
};

#endif // CORONARYPROFILE_H

