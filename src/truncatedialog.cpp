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

#include "truncatedialog.h"
#include "ui_truncatedialog.h"


TruncateDialog::TruncateDialog(QWidget* parent)
	:	QDialog(parent),
		ui(new Ui::TruncateDialog)
{
	ui->setupUi(this);

	connect(ui->sb_start, SIGNAL(valueChanged(double)), this, SLOT(startChanged(double)));
    connect(ui->sb_end, SIGNAL(valueChanged(double)),   this, SLOT(endChanged(double)));
}

TruncateDialog::~TruncateDialog()
{
	delete ui;
}

void TruncateDialog::setRange(double start, double end)
{
	ui->sb_start->setMinimum(start);
	ui->sb_start->setMaximum(end);
	ui->sb_start->setValue(start);
	
	ui->sb_end->setMinimum(start);
	ui->sb_end->setMaximum(end);
	ui->sb_end->setValue(end);
}

void TruncateDialog::startChanged(double val)
{
	ui->sb_end->setMinimum(val);
}

void TruncateDialog::endChanged(double val)
{
	ui->sb_start->setMaximum(val);
}

double TruncateDialog::getStart() const
{
	return ui->sb_start->value();
}

double TruncateDialog::getEnd() const
{
	return ui->sb_end->value();
}

double TruncateDialog::getStep() const
{
	return ui->sb_start->singleStep();
}

