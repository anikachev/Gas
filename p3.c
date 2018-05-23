#include <GL/glut.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#define N 2000 // число молекул

GLint Width = 700, Height = 620; //размеры окна

const int Pointsize = 3; // размер молекулы и её радиуса взаимодействия
const int Diametr = 3;
const float Pi = 3.1415926535;
const float MaxSpeedProj = 8; // максимальное значение каждой из проекций скоростей
float MaxSpeed;
int GraphConst = 0;
int PointConst = 1;
int grav = -1;
int Sten = -1;
unsigned long dat[300] = {0};
unsigned long tdat[100] = {0};
int MaxwellGraph[200] = {0};
int MaxwellCycle = 0;
int MaxwellData[200] = {0};


struct Particle
{
	float x[2];
	float vx[2];
	int mass;
	int psize;
};
struct Particle p[N]; // массив молекул

float speed(struct Particle p)
{
	float s = sqrt(p.vx[0]*p.vx[0] + p.vx[1]*p.vx[1]);
	return s; // измерение модуля скорости молекулы
}

float dist(struct Particle p1, struct Particle p2)
{
	float d = sqrt((p1.x[0] - p2.x[0])*(p1.x[0] - p2.x[0]) + (p1.x[1] - p2.x[1])*(p1.x[1] - p2.x[1]));
	return d;
}

void DrawPoint(struct Particle p) // рисует молекулу
{
	float s = speed(p);
	glColor3f(1.0*(s/MaxSpeed), 0, 1.0*(1-s/MaxSpeed)); // задаёт цвет молекулы в зависимости от модуля скорости

	glBegin(GL_QUADS);
		glVertex2f(p.x[0] - p.psize/2, p.x[1] - p.psize/2);
		glVertex2f(p.x[0] - p.psize/2, p.x[1] + p.psize/2);
		glVertex2f(p.x[0] + p.psize/2, p.x[1] + p.psize/2);
		glVertex2f(p.x[0] + p.psize/2, p.x[1] - p.psize/2);
	glEnd();
}

void DrawGraph0(struct Particle *p) // рисует график температуры от координаты
{
	int i, k, j;
	float v, d;
	glColor3f(0, 0, 0);

	d = Width/100.0;
	for(i = 1; i <= 100; ++i) // Находит средн. скорость на каждом из 100 участков
	{
		k = 0;
		v = 0;
		for(j = 0; j < N; ++j)
		{
			if(p[j].x[0] <= d*i && p[j].x[0] >= d*(i-1))
			{
				k++;
				v += speed(p[j]);
			}
		}
        if(k != 0)
			v = v/k;
		glBegin(GL_LINES);
			glVertex2f(d*(i-1), (d*v + 1));
			glVertex2f(d*i, (d*v + 1));
		glEnd();
		tdat[i] += v; // записывает в файл скорость на i-м участке
	}
}

void DrawGraph1(struct Particle *p) // график распределения Максвелла
{
	double d;
	int i, j, k = 0;
	d = 2*MaxSpeed/200.0; // делим возможные скорости(от 0 до 2*MaxSpeed) на 200 частей
	glColor3f(0, 0, 0);
	
	MaxwellCycle++;  
	if(MaxwellCycle >= 5)// усредняет результаты для плавного графика
		{
			for(i = 1; i <= 200; i++)
			{
				MaxwellData[i] = MaxwellGraph[i]/5;
				MaxwellGraph[i] = 0;
				//printf("%d\n", MaxwellGraph[50]);
			}
			MaxwellCycle = 0;
		}
	//printf("%d\n%d\n\n", MaxwellGraph[50], MaxwellData[50]);
	for(i = 1; i <= 200; ++i) 
	{
		k = 0;
		for(j = 0; j < N; ++j)// для каждого участка скорости находим число молекул
		{
			if(speed(p[j]) > (i-1)*d && speed(p[j]) < i*d)
				k++;
		}
		dat[i] += k; // записываем в файл		
		//printf("%d\n", k);
		MaxwellGraph[i] += k;
		glBegin(GL_LINES);
			glVertex2f(Width/200.0*(i-1), 30*MaxwellData[i]*Height/(double)N);
			glVertex2f(Width/200.0*i, 30*MaxwellData[i]*Height/(double)N);
		glEnd();
	}
}
float stenka(int k) // генерирует рандомные скорости от 2k/3 до 4k/3;
{
	float h;
	h = rand()/(float)RAND_MAX*2*k/3 + 2*k/3;
	return h;
}

void Generate(struct Particle *_p, int n, int mass, int size, float _x, float _y, int lengx, int lengy, int deg, int speed)
{
	long i;
	for(i = 0; i < n; ++i)
		{
			_p[i].x[0] = rand() % lengx + _x;
			_p[i].x[1] = rand() % lengy + _y + Height/5;
			
			_p[i].vx[0] = stenka(speed*cos((double)deg/180*Pi));
			_p[i].vx[1] = stenka(speed*sin((double)deg/180*Pi));
			
			_p[i].psize = size;
			_p[i].mass = mass;
		}
}

void Gravity(struct Particle *_p)
{
	int i;
	for(i = 0; i < N; i++)
	{
			_p[i].vx[1] -= 1;
	}
}
void InputFunc()
{
		int n, mass, size, dx, dy, deg, speed;
		int k = 0;
		float x, y;
		FILE * in = fopen("input.txt", "r");
		while(!feof(in))
		{
				if(fscanf(in, "%d %d %d %f %f %d %d %d %d", &n, &mass, &size, &x, &y, &dx, &dy, &deg, &speed))
				{
					k += n;
					Generate(p + k, n, mass, size, x, y, dx, dy, deg, speed);
				}
		}
		fclose(in);
}
					
void Display(void)
{
	int i;
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);

	if(PointConst == 1) // управление отображением молекул
		for(i = 0; i < N; i++)
			DrawPoint(p[i]);
			
	glColor3f(1, 1, 1); // генерирует окно для графика
	glBegin(GL_QUADS);
		glVertex2f(0, 0);
		glVertex2f(0, Height/5);
		glVertex2f(Width, Height/5);
		glVertex2f(Width, 0);
	glEnd();

	if(GraphConst == 0) // переключение графиков
		DrawGraph0(p);

	if(GraphConst == 1)
		DrawGraph1(p);

	glFinish();
}

void Reshape(GLint w, GLint h)
{
	Width = w;
	Height = h;

	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, w, 0, h, -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void tick()
{
	int i, j;
	float temp, mom[2], nmom[2], tan[2], ntan[2], v1tan, v2tan, v1nor, v2nor;
	if(Sten == -1) // горячая стенка отключена
		for(i = 0; i < N; ++i) // проверка на столкновения со стеной
		{
			if(p[i].x[0] > Width - p[i].psize/2)
				if(p[i].vx[0] > 0)
					p[i].vx[0] *= -1;
			if(p[i].x[0] < p[i].psize/2)
				if(p[i].vx[0] < 0)
					p[i].vx[0] *= -1;
			if(p[i].x[1] > Height - p[i].psize/2)
				if(p[i].vx[1] > 0)
					p[i].vx[1] *= -1;
			if(p[i].x[1] < (p[i].psize + Height/5))
				if(p[i].vx[1] < 0)
					p[i].vx[1] *= -1;
		}
	else // если горячая стенка включена
		for(i = 0; i < N; ++i)
		{
			if(p[i].x[0] > Width - p[i].psize/2)
				if(p[i].vx[0] > 0)
					p[i].vx[0] = stenka(-2);
			if(p[i].x[0] < p[i].psize/2)
				if(p[i].vx[0] < 0)
					p[i].vx[0] = stenka(10);
			if(p[i].x[1] > Height - p[i].psize/2)
				if(p[i].vx[1] > 0)
					p[i].vx[1] *= -1;
			if(p[i].x[1] < (p[i].psize/2 + Height/5))
				if(p[i].vx[1] < 0)
					p[i].vx[1] *= -1;
		}

	for(i = 0; i < N - 1; ++i) // обработка столкновений частиц
	{
		for(j = i + 1; j < N; ++j)
		{
			if(dist(p[i], p[j]) < (p[i].psize + p[j].psize)/2.0)
			{
				mom[0] = p[j].x[0] - p[i].x[0];
				mom[1] = p[j].x[1] - p[i].x[1];
				tan[0] = mom[1];
				tan[1] = -mom[0];

				ntan[0] = tan[0]/dist(p[i], p[j]);
				ntan[1] = tan[1]/dist(p[i], p[j]);

				nmom[0] = mom[0]/dist(p[i], p[j]);
				nmom[1] = mom[1]/dist(p[i], p[j]);

				//printf("%f\t,%f\n\n%f\t,%f\n", ntan[0], ntan[1], nmom[0], nmom[1]);


				v1tan = p[i].vx[0]*ntan[0] + p[i].vx[1]*ntan[1];
				v2tan = p[j].vx[0]*ntan[0] + p[j].vx[1]*ntan[1];

				v1nor = p[i].vx[0]*nmom[0] + p[i].vx[1]*nmom[1];
				v2nor = p[j].vx[0]*nmom[0] + p[j].vx[1]*nmom[1];

				temp = (2*p[i].mass*v1nor + v2nor*(p[j].mass - p[i].mass))/(p[i].mass + p[j].mass);
				v1nor = (2*p[j].mass*v2nor + v1nor*(p[i].mass - p[j].mass))/(p[i].mass + p[j].mass);
				v2nor = temp;

				p[i].vx[0] = ntan[0]*v1tan + nmom[0]*v1nor;
				p[i].vx[1] = ntan[1]*v1tan + nmom[1]*v1nor;

				p[j].vx[0] = ntan[0]*v2tan + nmom[0]*v2nor;
				p[j].vx[1] = ntan[1]*v2tan + nmom[1]*v2nor;
			}
		}

	}

	for(i = 0; i < N; ++i) // движение частиц
	{
		p[i].x[0] += p[i].vx[0];
		p[i].x[1] += p[i].vx[1];
	}
	
	if(grav == 1)
		Gravity(p);
}

void timf(int a)
{
	tick();
	glutPostRedisplay();
	glutTimerFunc(10, timf, 0);
	//printf("%lu\n%lu\n%lu\n%lu\n\n", dat[0], dat[1], dat[2], dat[3]);
}

void KBFUNK(unsigned char g, int x, int y)
{
	int i;
	float d;
	FILE * out = fopen("out.txt", "a");
	FILE * out2 = fopen("out2.txt", "w");
	switch(g)
	{
		case '0':
			GraphConst = 0;
			break;
		case '1':
			GraphConst = 1;
			break;
		case 'p':
			PointConst *= -1;
			break;
		case 'd':
			d = 2*MaxSpeed/200.0;
			for(i = 0; i < 200; ++i)
				fprintf(out, "%f\t%f\n", i*d, dat[i]/((double)N));
			for(i = 0; i < 100; ++i)
				fprintf(out2, "%u\t%lu\n", i, tdat[i]);

			break;
		case 'r':
			//Generate(p, N/2, 1, 3, 0, 0, Width, Height);
			//Generate(p + N/2, N/2, 2, 4,  0, 0, Width, Height);
			break;
		case 'u':
			Generate(p, 1000, 1, 5, 0, 0, Width, Height, 0, 8);
			break;
		case 'g':
			grav *= -1;
			break;
		case 's':
			Sten *= -1;
			break;
		case 'c':
			for(i = 0; i < 300; i++)
				dat[i] = 0;
			for(i = 0; i < 100; i++)
				tdat[i] = 0;
			break;
		case 'f':
			InputFunc();
			break;
		case 'o':
			for(i = 0; i < N; i++)
					p[i].x[0] = p[i].x[1] = p[i].vx[0] = p[i].vx[1] = 0;
			break;
	}
	fclose(out);
	fclose(out2);
}

int main(int argc, char *argv[])
{
	int i;
	MaxSpeed = MaxSpeedProj*sqrt(2);
	srand(time(NULL));
	//Generate(p, N);

	glutInit(&argc, argv);
   	glutInitDisplayMode(GLUT_RGB);
   	glutInitWindowSize(Width, Height);
   	glutCreateWindow("Single Particle");

   	glutDisplayFunc(Display);
   	glutReshapeFunc(Reshape);
	glutTimerFunc(40, timf, 0);
	glutKeyboardFunc(KBFUNK);

   	glutMainLoop();

}
