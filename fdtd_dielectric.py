# -*- coding: utf-8 -*-
'''
Моделирование отражения гармонического сигнала от слоя диэлектрика
'''

import math

import numpy

import tools
import sources
import boundary
from objects import LayerContinuous, LayerDiscrete, Probe


class Sampler:
    def __init__(self, discrete: float):
        self.discrete = discrete

    def sample(self, x: float) -> int:
        return math.floor(x / self.discrete + 0.5)


def sampleLayer(layer_cont, sampler):
    start_discrete = sampler.sample(layer_cont.xmin)
    end_discrete = (sampler.sample(layer_cont.xmax)
                    if layer_cont.xmax is not None
                    else None)
    return LayerDiscrete(start_discrete, end_discrete,
                         layer_cont.eps, layer_cont.mu)


def fillMedium(layer: LayerDiscrete, eps, mu):
    if layer.xmax is not None:
        eps[layer.xmin: layer.xmax] = layer.eps
        mu[layer.xmin: layer.xmax] = layer.mu
    else:
        eps[layer.xmin:] = layer.eps
        mu[layer.xmin:] = layer.mu


if __name__ == '__main__':
    # Используемые константы
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Скорость света в вакууме
    c = 299792458.0

    # Электрическая постоянная
    eps0 = 8.854187817e-12

    # Параметры моделирования
    # Частота сигнала, Гц
    f_Hz = 0.2e9

    # Дискрет по пространству в м
    dx = 0.07495

    wavelength = c / f_Hz
    Nl = wavelength / dx

    # Число Куранта
    Sc = 1.0

    # Размер области моделирования в м
    maxSize_m = 14.99

    # Время расчета в секундах
    maxTime_s = 7.999e-8

    # Положение источника в м
    sourcePos_m = 3.7475

    # Координаты датчиков для регистрации поля в м
    probesPos_m = [1.87375, 3.7475, 11.2425]

    # Параметры слоев
    layers_cont = [LayerContinuous(7.495, eps=2.56)]

    # Скорость обновления графика поля
    speed_refresh = 30

    # Переход к дискретным отсчетам
    # Дискрет по времени
    dt = dx * Sc / c

    sampler_x = Sampler(dx)
    sampler_t = Sampler(dt)

    # Время расчета в отсчетах
    maxTime = sampler_t.sample(maxTime_s)

    # Размер области моделирования в отсчетах
    maxSize = sampler_x.sample(maxSize_m)

    # Положение источника в отсчетах
    sourcePos = sampler_x.sample(sourcePos_m)

    layers = [sampleLayer(layer, sampler_x) for layer in layers_cont]

    # Датчики для регистрации поля
    probesPos = [sampler_x.sample(pos) for pos in probesPos_m]
    probes = [Probe(pos, maxTime) for pos in probesPos]

    # Вывод параметров моделирования
    print('Число Куранта: {}'.format(Sc))
    print('Размер области моделирования: {} м'.format(maxSize_m))
    print('Время расчета: {} нс'.format(maxTime_s * 1e9))
    print('Координата источника: {} м'.format(sourcePos_m))
    print('Частота сигнала: {} ГГц'.format(f_Hz * 1e-9))
    print('Длина волны: {} м'.format(wavelength))
    print('Количество отсчетов на длину волны (Nl): {}'.format(Nl))
    probes_m_str = ', '.join(['{:.6f}'.format(pos) for pos in probesPos_m])
    print('Дискрет по пространству: {} м'.format(dx))
    print('Дискрет по времени: {} нс'.format(dt * 1e9))
    print('Координаты пробников [м]: {}'.format(probes_m_str))
    print()
    print('Размер области моделирования: {} отсч.'.format(maxSize))
    print('Время расчета: {} отсч.'.format(maxTime))
    print('Координата источника: {} отсч.'.format(sourcePos))
    probes_str = ', '.join(['{}'.format(pos) for pos in probesPos])
    print('Координаты пробников [отсч.]: {}'.format(probes_str))

    # Параметры среды
    # Диэлектрическая проницаемость
    eps = numpy.ones(maxSize)

    # Магнитная проницаемость
    mu = numpy.ones(maxSize - 1)

    for layer in layers:
        fillMedium(layer, eps, mu)

    # Источник
    amp = 1.0
    amp_E = amp * Sc / numpy.sqrt(eps[sourcePos] * mu[sourcePos])
    amp_H = amp * Sc / (W0 * mu[sourcePos])

    source_E = sources.make_harmonic(amp_E, f_Hz, Sc, dt)
    source_H = sources.make_harmonic(amp_H, f_Hz, Sc, dt)

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize - 1)

    # Создание экземпляров классов граничных условий
    boundary_left = boundary.ABCSecondLeft(eps[0], mu[0], Sc)
    boundary_right = boundary.ABCSecondRight(eps[-1], mu[-1], Sc)

    # Параметры отображения поля E
    display_field = Ez
    display_ylabel = 'Ez, В/м'
    display_ymin = -2.1
    display_ymax = 2.1

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
    display = tools.AnimateFieldDisplay(dx, dt,
                                        maxSize,
                                        display_ymin, display_ymax,
                                        display_ylabel,
                                        title='fdtd_dielectric')

    display.activate()
    display.drawSources([sourcePos])
    display.drawProbes(probesPos)
    for layer in layers:
        display.drawBoundary(layer.xmin)
        if layer.xmax is not None:
            display.drawBoundary(layer.xmax)

    for t in range(1, maxTime):
        # Расчет компоненты поля H
        Hy = Hy + (Ez[1:] - Ez[:-1]) * Sc / (W0 * mu)

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[sourcePos - 1] -= source_H.getField(t)

        # Расчет компоненты поля E
        Ez[1:-1] = Ez[1: -1] + (Hy[1:] - Hy[: -1]) * Sc * W0 / eps[1: -1]

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[sourcePos] += source_E.getField(t + 0.5 - (-0.5))

        boundary_left.updateField(Ez, Hy)
        boundary_right.updateField(Ez, Hy)

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if t % speed_refresh == 0:
            display.updateData(display_field, t)

    display.stop()

    # Отображение сигнала, сохраненного в пробнике
    tools.showProbeSignals(probes, dx, dt, -2.1, 2.1)
