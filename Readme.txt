Readme
Данная программа позволяет проводить моделирование физических процессов с идеальным газом.
Список горячих клавиш:
P - включает/отключает отображение молекул на экране. (может помочь производительности при большом числе молекул) 

S - включает/отключает режим, при котором левая стенка сосуда нагревается, а правая охлаждается (может быть изучено распределение температуры по ширине сосуда)

Функции “0”, “1”, ... (to be continued) управляют графиком под окном вывода:
1 - построение зависимости вероятности обнаружения молекулы от модуля скорости( распределение Максвелла)
0 - график значения средней температуры молекул на данном участке ширины сосуда 

Для вывода более полной информации о распределениях в текстовый файл нажмите клавишу
D - создает файлы οut.txt, out2.txt… 
Описание файлов:
1. out.txt - 
2. out2.txt -

Для включения притяжения нажмите клавишу g(тестовая функция!!! последствия необратимы!!!)

Для передачи параметров частиц создайте файл input.txt 
Каждая строчка задаёт семейство частиц с параметрами через пробел 
: кол-во, масса, диаметр, координаты прямоугольника x, y, ширина и высота dx, dy, угол скоростей, средн. скорость.
Для запуска файла нажмите F

Описание возможностей исходного кода:

Функция Generate(...) позволяет сгенерировать N частиц(помещаются в адрес p) массы m, диаметра d, со средней скоростью ν(распространение во все стороны равновероятно) в прямоугольнике, построенном на точках (x, y),  (x + Δx, y + Δy).


