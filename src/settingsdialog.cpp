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

#include "settingsdialog.h"
#include "ui_settingsdialog.h"

#include "settings.h"

#include <QFileDialog>


SettingsDialog::SettingsDialog(QWidget* parent):
	QDialog(parent),
	ui(new Ui::SettingsDialog)
{
	ui->setupUi(this);

	connect(ui->b_browse, SIGNAL(clicked()), this, SLOT(browse_clicked()));

	ui->tabWidget->setCurrentIndex(0);
}

SettingsDialog::~SettingsDialog()
{
	delete ui;
}

void SettingsDialog::browse_clicked()
{
	QString dir = QFileDialog::getExistingDirectory(this, "Select Directory", ui->defaultDir->text());
	if (!dir.isEmpty())
		ui->defaultDir->setText(dir);
}

void SettingsDialog::setupFromConfig()
{
	Settings& s = Settings::getSettings();

	setNumProfiles(s.getNumNormalsAtPoint());
	setSizeProfile(s.getSizeNormal());
	setMinRelMagnitude(s.getMinRelMagniude());
	setPlateauThreshold(s.getPlateauThreshold());
	setMinVesselDiameter(s.getMinVesselDiameter());
	setSplineSampling(s.getSharpnessEvalInterval());
	setCheckPoints(s.getDoPointCheck());
	setCheckNormals(s.getDoNormalCheck());
	setUseEllipse(s.getUseEllipse());
	setShowMIP(s.getShowMIP());
	setDefaultDir(s.getDefaultDir());
	setRememberLastDir(s.getRememberLastDir());
}

void SettingsDialog::updateConfig() const
{
	Settings& s = Settings::getSettings();

	s.setNumNormalsAtPoint(getNumProfiles());
	s.setSizeNormal(getSizeProfile());
	s.setMinRelMagnitude(getMinRelMagnitude());
	s.setPlateauThreshold(getPlateauThreshold());
	s.setMinVesselDiameter(getMinVesselDiameter());
	s.setSharpnessEvalInterval(getSplineSampling());
	s.setDoPointCheck(getCheckPoints());
	s.setDoNormalCheck(getCheckNormals());
	s.setUseEllipse(getUseEllipse());
	s.setShowMIP(getShowMIP());
	s.setDefaultDir(getDefaultDir());
	s.setRememberLastDir(getRememberLastDir());
}

void SettingsDialog::setNumProfiles(size_t val)
{
	ui->sb_numProfiles->setValue(val);
}

void SettingsDialog::setSizeProfile(size_t val)
{
	ui->sb_sizeProfiles->setValue(val);
}

void SettingsDialog::setMinRelMagnitude(double val)
{
	ui->sb_minRelMagnitude->setValue(val);
}

void SettingsDialog::setPlateauThreshold(double val)
{
	ui->sb_plateauThreshold->setValue(val);
}

void SettingsDialog::setMinVesselDiameter(double val)
{
	ui->sb_minDiameter->setValue(val);
}

void SettingsDialog::setSplineSampling(double val)
{
	ui->sb_splineSampling->setValue(val);
}

void SettingsDialog::setCheckPoints(bool val)
{
	ui->cb_checkPoints->setChecked(val);
}

void SettingsDialog::setCheckNormals(bool val)
{
	ui->cb_checkNormals->setChecked(val);
}

void SettingsDialog::setUseEllipse(bool val)
{
	ui->cb_useEllipse->setChecked(val);
}

void SettingsDialog::setShowMIP(bool val)
{
	ui->cb_mip->setChecked(val);
}

void SettingsDialog::setDefaultDir(const QString& val)
{
	ui->defaultDir->setText(val);
}

void SettingsDialog::setRememberLastDir(bool val)
{
	ui->cb_rememberLastDir->setChecked(val);
}

size_t SettingsDialog::getNumProfiles() const
{
	return ui->sb_numProfiles->value();
}

size_t SettingsDialog::getSizeProfile() const
{
	return ui->sb_sizeProfiles->value();
}

double SettingsDialog::getMinRelMagnitude() const
{
	return ui->sb_minRelMagnitude->value();
}

double SettingsDialog::getPlateauThreshold() const
{
	return ui->sb_plateauThreshold->value();
}

double SettingsDialog::getMinVesselDiameter() const
{
	return ui->sb_minDiameter->value();
}

double SettingsDialog::getSplineSampling() const
{
	return ui->sb_splineSampling->value();
}

bool SettingsDialog::getCheckPoints() const
{
	return ui->cb_checkPoints->isChecked();
}

bool SettingsDialog::getCheckNormals() const
{
	return ui->cb_checkNormals->isChecked();
}

bool SettingsDialog::getUseEllipse() const
{
	return ui->cb_useEllipse->isChecked();
}

bool SettingsDialog::getShowMIP() const
{
	return ui->cb_mip->isChecked();
}

QString SettingsDialog::getDefaultDir() const
{
	return ui->defaultDir->text();
}

bool SettingsDialog::getRememberLastDir() const
{
	return ui->cb_rememberLastDir->isChecked();
}

