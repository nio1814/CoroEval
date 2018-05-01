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

#ifndef TRUNCATEDIALOG_H
#define TRUNCATEDIALOG_H

#include <QDialog>


namespace Ui {
	class TruncateDialog;
}

class TruncateDialog : public QDialog
{
	Q_OBJECT

public:
	explicit TruncateDialog(QWidget* parent = 0);
	~TruncateDialog();

	void setRange(double start, double end);
	double getStart() const;
	double getEnd() const;
	double getStep() const;

private:
	Ui::TruncateDialog* ui;

private slots:
	void startChanged(double val);
    void endChanged(double val);
};

#endif // TRUNCATEDIALOG_H

