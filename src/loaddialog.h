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

#ifndef LOADDIALOG_H
#define LOADDIALOG_H

#include <QDialog>


// Forward declarations
class Settings;
class Data;

namespace Ui {
class LoadDialog;
}

namespace gdcm {
class PixelFormat;
class Image;
}


class LoadDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit LoadDialog(QWidget* parent = 0);
    ~LoadDialog();

	void setLastLoadDir(const QString& dir);
	QString getLastLoadDir() const;

protected slots:
    void browse_clicked();
    void loadFile();
	void cb_changed(int state);

signals:
	void clearData() const;
    void sendData(Data*) const;
    void sendFileName(QString) const;

private:
	Settings& m_settings;
    Ui::LoadDialog* ui;

	QString m_lastLoadDir;

	void loadRaw();
	void loadCavarev();
	void loadDicomDir();
	void loadDicomZIP();
	bool loadDicomSlice(const gdcm::PixelFormat& pixeltype, const gdcm::Image& image, size_t sliceOffset, size_t index, float* data);
};

#endif // LOADDIALOG_H

