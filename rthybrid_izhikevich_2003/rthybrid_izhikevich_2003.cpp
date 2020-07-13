/*
 * Copyright (C) 2011 Georgia Institute of Technology, University of Utah,
 * Weill Cornell Medical College
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * This is a template implementation file for a user module derived from
 * DefaultGUIModel with a custom GUI.
 */

#include "rthybrid_izhikevich_2003.h"
#include <iostream>
#include <main_window.h>



#define NM_IZHIKEVICH_2003_DT 0
#define NM_IZHIKEVICH_2003_I 1
#define NM_IZHIKEVICH_2003_SYN 2
#define NM_IZHIKEVICH_2003_A 3
#define NM_IZHIKEVICH_2003_B 4
#define NM_IZHIKEVICH_2003_C 5
#define NM_IZHIKEVICH_2003_D 6

//Vars
#define NM_IZHIKEVICH_2003_V 0
#define NM_IZHIKEVICH_2003_U 1

extern "C" Plugin::Object*
createRTXIPlugin(void)
{
  return new RTHybridIzhikevich2003();
}

static DefaultGUIModel::variable_t vars[] = {
    {"Burst duration (s)", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    {"I", "",DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

    { "a", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "b", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "c", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "d", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},


    { "v0", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
    { "u0", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},

    { "Vm (v)", "Membrane potential (in V)", DefaultGUIModel::OUTPUT, },
    { "Vm (mV)", "Membrane potential (in mV)", DefaultGUIModel::OUTPUT, },

    { "Isyn (nA)", "Synaptic input current (in nA)", DefaultGUIModel::INPUT,},
    { "Burst duration (s)", "", DefaultGUIModel::INPUT,},

    { "v", "", DefaultGUIModel::STATE,},
    { "s_points", "", DefaultGUIModel::STATE,},
    { "dt", "", DefaultGUIModel::STATE,},
    { "syn", "", DefaultGUIModel::STATE,},
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

RTHybridIzhikevich2003::RTHybridIzhikevich2003(void)
  : DefaultGUIModel("RTHybrid Izhikevich (2003) neuron model", ::vars, ::num_vars)
{
  setWhatsThis("<p><b>RTHybrid Izhikevich (2003) neuron model</b></p>");
  DefaultGUIModel::createGUI(vars,
                             num_vars); // this is required to create the GUI
  //customizeGUI();
  initParameters();
  update(INIT); // this is optional, you may place initialization code directly
                // into the constructor
  refresh();    // this is required to update the GUI with parameter and state
                // values
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

RTHybridIzhikevich2003::~RTHybridIzhikevich2003(void)
{
}


void RTHybridIzhikevich2003::runge_kutta_65 (void (*f) (double *, double *, double *, double), int dim, double dt, double * vars, double * params, double aux) {
    double apoyo[dim], retorno[dim];
    double k[6][dim];
    int j;

    (*f)(vars, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[0][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.2;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[1][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.075 + k[1][j] * 0.225;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[2][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.3 - k[1][j] * 0.9 + k[2][j] * 1.2;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[3][j] = dt * retorno[j];
        apoyo[j] = vars[j] + k[0][j] * 0.075 + k[1][j] * 0.675 - k[2][j] * 0.6 + k[3][j] * 0.75;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[4][j] = dt * retorno[j];
        apoyo[j] = vars[j]
                + k[0][j] * 0.660493827160493
                + k[1][j] * 2.5
                - k[2][j] * 5.185185185185185
                + k[3][j] * 3.888888888888889
                - k[4][j] * 0.864197530864197;
    }

    (*f)(apoyo, retorno, params, aux);
    for(j = 0; j < dim; ++j) {
        k[5][j] = dt * retorno[j];
        apoyo[j] = vars[j]
                + k[0][j]*0.1049382716049382
                + k[2][j]*0.3703703703703703
                + k[3][j]*0.2777777777777777
                + k[4][j]*0.2469135802469135;
    }


    for(j = 0; j < dim; ++j) {
        vars[j] += k[0][j]*0.098765432098765+
                   k[2][j]*0.396825396825396+
                   k[3][j]*0.231481481481481+
                   k[4][j]*0.308641975308641-
                   k[5][j]*0.035714285714285;
    }

    return;
}


void RTHybridIzhikevich2003::select_dt_neuron_model (double * dts, double * pts, unsigned int length, double pts_live, double * dt, double * pts_burst) {
    double aux = pts_live;
    double factor = 1;
    double intpart, fractpart;
    int flag = 0;
    int i;

    *dt = -1;
    *pts_burst = -1;

    while (aux < pts[0]) {
        aux = pts_live * factor;
        factor += 1;

        for (i = length - 1; i >= 0; i--) {
          if (pts[i] > aux) {
            *dt = dts[i];
              *pts_burst = pts[i];

              fractpart = modf(*pts_burst / pts_live, &intpart);

              if (fractpart <= 0.1*intpart) flag = 1;

              break;
          }
            
        }

        if (flag == 1) break;
    }

    if (flag == 0) {
        for (i = length - 1; i >= 0; i--) {
          if (pts[i] > aux) {
            *dt = dts[i];
              *pts_burst = pts[i];

              break;  
          }
        }
    }

    return;
}


double RTHybridIzhikevich2003::set_pts_burst (double sec_per_burst) {
  int length = 0;
  int method = 3;
  double pts_match = sec_per_burst * freq;
  double pts_burst, dt;


  length =  83.000000 - 19;
  double dts[] = {0.000010, 0.000200, 0.000300, 0.000400, 0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000, 0.001100, 0.001200, 0.001300, 0.001400, 0.001500, 0.001600, 0.001700, 0.001800, 0.001900, 0.002000, 0.002100, 0.002200, 0.002300, 0.002400, 0.002600, 0.002800, 0.003000, 0.003200, 0.003400, 0.003700, 0.004000, 0.004400, 0.004800, 0.005300, 0.005900, 0.006000, 0.006100, 0.006200, 0.006300, 0.006400, 0.006500, 0.006600, 0.006700, 0.006800, 0.006900, 0.007000, 0.007100, 0.007200, 0.007300, 0.007400, 0.007500, 0.007600, 0.007800, 0.008000, 0.008200, 0.008400, 0.008600, 0.008800, 0.009000, 0.009200, 0.009400, 0.009600, 0.009800, 0.010000/*, 0.010200, 0.010400, 0.010600, 0.010800, 0.011100, 0.011400, 0.011700, 0.012000, 0.012300, 0.012600, 0.012900, 0.092600, 0.096900, 0.097000, 0.097200, 0.097400, 0.097600, 0.098200, 0.100000*/};
  double pts[] = {383275.000000, 296562.000000, 197709.000000, 148282.800000, 118627.285714, 98856.666667, 84736.000000, 74143.500000, 65907.071429, 59317.933333, 53924.941176, 49430.421053, 45630.350000, 42371.000000, 39544.625000, 37075.080000, 34893.666667, 32954.172414, 31222.645161, 29662.000000, 28247.000000, 26965.555556, 25792.945946, 24716.153846, 22817.023810, 21187.673913, 19773.979592, 18540.673077, 17447.125000, 16033.983607, 14832.015152, 13483.000000, 12361.000000, 11195.000000, 10059.500000, 9888.030000, 9726.009901, 9572.495146, 9418.990476, 9274.009434, 9129.962963, 8992.000000, 8857.477477, 8726.000000, 8601.669565, 8478.258621, 8359.000000, 8240.000000, 8129.000000, 8017.032520, 7912.504000, 7810.000000, 7609.676923, 7419.015038, 7235.000000, 7064.992857, 6899.000000, 6743.496599, 6594.006667, 6450.480519, 6313.496815, 6183.000000, 6056.500000, 5936.982036/*, 5820.476471, 5708.011494, 5598.627119, 5495.500000, 5349.010811, 5208.005236, 5074.000000, 4946.000000, 4827.004854, 4709.758294, 4603.500000, 4497.207207, 2146.965591, 1043.224660, 842.464979, 783.403137, 674.209177, 610.585575, 599.731092*/};
  select_dt_neuron_model(dts, pts, length, pts_match, &(params_model[NM_IZHIKEVICH_2003_DT]), &(pts_burst));

  printf("%f %f\n", pts_burst, params_model[NM_IZHIKEVICH_2003_DT]);

  return pts_burst;
}


void RTHybridIzhikevich2003::nm_izhikevich_2003_f (double * vars, double * ret, double * params, double syn) {
    ret[NM_IZHIKEVICH_2003_V] = 0.04 * vars[NM_IZHIKEVICH_2003_V]*vars[NM_IZHIKEVICH_2003_V] + 5*vars[NM_IZHIKEVICH_2003_V] + 140 - vars[NM_IZHIKEVICH_2003_U] + params[NM_IZHIKEVICH_2003_I] - syn;
    ret[NM_IZHIKEVICH_2003_U] = params[NM_IZHIKEVICH_2003_A] * (params[NM_IZHIKEVICH_2003_B] * vars[NM_IZHIKEVICH_2003_V] - vars[NM_IZHIKEVICH_2003_U]);

    return;
}


void
RTHybridIzhikevich2003::execute(void)
{
    int i;
    for (i = 0; i < s_points; i++) {
      runge_kutta_65(nm_izhikevich_2003_f, 2, 0.001, vars_model, params_model, input(0));
    }

    if (vars_model[NM_IZHIKEVICH_2003_V] >= 30) {
        vars_model[NM_IZHIKEVICH_2003_V] = params_model[NM_IZHIKEVICH_2003_C];
        vars_model[NM_IZHIKEVICH_2003_U] += params_model[NM_IZHIKEVICH_2003_D];
    }

    output(0) = vars_model[NM_IZHIKEVICH_2003_V] / 1000.0;
    output(1) = vars_model[NM_IZHIKEVICH_2003_V];

    return;
}

void
RTHybridIzhikevich2003::initParameters(void)
{
    burst_duration = 1.0;
    freq = 1.0 / (period * 1e-3);
    s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq)); 
    if (s_points == 0) s_points = 1;

    params_model[NM_IZHIKEVICH_2003_A] = 0.02;
    params_model[NM_IZHIKEVICH_2003_B] = 0.2;
    params_model[NM_IZHIKEVICH_2003_C] = -50.0;
    params_model[NM_IZHIKEVICH_2003_D] = 2.0;
    params_model[NM_IZHIKEVICH_2003_I] = 10.0;
    vars_model[NM_IZHIKEVICH_2003_V] = 30.2403;
    vars_model[NM_IZHIKEVICH_2003_U] = -5.54459;
}

void
RTHybridIzhikevich2003::update(DefaultGUIModel::update_flags_t flag)
{
  switch (flag) {
    case INIT:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      freq = 1.0 / (period * 1e-3);

      setParameter("Burst duration (s)", burst_duration);
      setParameter("I", params_model[NM_IZHIKEVICH_2003_I]);

      setParameter("a", params_model[NM_IZHIKEVICH_2003_A]);
      setParameter("b", params_model[NM_IZHIKEVICH_2003_B]);
      setParameter("c", params_model[NM_IZHIKEVICH_2003_C]);
      setParameter("d", params_model[NM_IZHIKEVICH_2003_D]);
      
      setParameter("v0", vars_model[NM_IZHIKEVICH_2003_V]);
      setParameter("u0", vars_model[NM_IZHIKEVICH_2003_U]);

      setState("v", vars_model[NM_IZHIKEVICH_2003_V]);
      setState("s_points", s_points);
      setState("dt", params_model[NM_IZHIKEVICH_2003_DT]);
      setState("syn", params_model[NM_IZHIKEVICH_2003_SYN]);
      break;

    case MODIFY:
      burst_duration = getParameter("Burst duration (s)").toDouble();
      freq = 1.0 / (period * 1e-3);
      if (burst_duration <= -1) burst_duration = input(1);
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq));
      if (s_points < 1) s_points = 1;

      params_model[NM_IZHIKEVICH_2003_A] = getParameter("a").toDouble();
      params_model[NM_IZHIKEVICH_2003_B] = getParameter("b").toDouble();
      params_model[NM_IZHIKEVICH_2003_C] = getParameter("c").toDouble();
      params_model[NM_IZHIKEVICH_2003_D] = getParameter("d").toDouble();
      params_model[NM_IZHIKEVICH_2003_I] = getParameter("I").toDouble();

      vars_model[NM_IZHIKEVICH_2003_V] = getParameter("v0").toDouble();
      vars_model[NM_IZHIKEVICH_2003_U] = getParameter("u0").toDouble();
      break;

    case UNPAUSE:
      freq = 1.0 / (period * 1e-3);
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq)); 
      if (s_points == 0) s_points = 1;

      break;

    case PAUSE:
      break;

    case PERIOD:
      period = RT::System::getInstance()->getPeriod() * 1e-6; // ms
      freq = 1.0 / (period * 1e-3);
      s_points = (int)(set_pts_burst(burst_duration) / (burst_duration * freq)); 
      if (s_points == 0) s_points = 1;


      break;

    default:
      break;
  }
}

void
RTHybridIzhikevich2003::customizeGUI(void)
{
  QGridLayout* customlayout = DefaultGUIModel::getLayout();

  QGroupBox* button_group = new QGroupBox;

  QPushButton* abutton = new QPushButton("Button A");
  QPushButton* bbutton = new QPushButton("Button B");
  QHBoxLayout* button_layout = new QHBoxLayout;
  button_group->setLayout(button_layout);
  button_layout->addWidget(abutton);
  button_layout->addWidget(bbutton);
  QObject::connect(abutton, SIGNAL(clicked()), this, SLOT(aBttn_event()));
  QObject::connect(bbutton, SIGNAL(clicked()), this, SLOT(bBttn_event()));

  customlayout->addWidget(button_group, 0, 0);
  setLayout(customlayout);
}

// functions designated as Qt slots are implemented as regular C++ functions
void
RTHybridIzhikevich2003::aBttn_event(void)
{
}

void
RTHybridIzhikevich2003::bBttn_event(void)
{
}
