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

#ifndef SETTINGSDIALOG_H
#define SETTINGSDIALOG_H

#include <QDialog>


namespace Ui {
	class SettingsDialog;
}


class SettingsDialog : public QDialog
{
	Q_OBJECT

public:
	explicit SettingsDialog(QWidget* parent = 0);
	~SettingsDialog();

	void setupFromConfig();
	void updateConfig() const;

	void setNumProfiles(size_t val);
	void setSizeProfile(size_t val);
	void setMinRelMagnitude(double val);
	void setPlateauThreshold(double val);
	void setMinVesselDiameter(double val);
	void setSplineSampling(double val);
	void setCheckPoints(bool val);
	void setCheckNormals(bool val);
	void setUseEllipse(bool val);
	void setShowMIP(bool val);
	void setDefaultDir(const QString& val);
	void setRememberLastDir(bool val);

	size_t getNumProfiles() const;
	size_t getSizeProfile() const;
	double getMinRelMagnitude() const;
	double getPlateauThreshold() const;
	double getMinVesselDiameter() const;
	double getSplineSampling() const;
	bool getCheckPoints() const;
	bool getCheckNormals() const;
	bool getUseEllipse() const;
	bool getShowMIP() const;
	QString getDefaultDir() const;
	bool getRememberLastDir() const;

protected slots:
    void browse_clicked();

private:
	Ui::SettingsDialog* ui;
};

#endif // SETTINGSDIALOG_H

