// Modifictaion of https://www.shadertoy.com/view/4sc3Wn

float pi = atan(1.0)*4.0;
float tau = atan(1.0)*8.0;

float scale = 1.0 / 9.0;

float epsilon = 1e-3;
float infinity = 1e6;

//Settings
//Uses cheaper arcs for common sweep angles (90 & 180 degrees).
#define USE_CHEAP_ARCS

#define TEXT_COLOR   vec3(0.00, 0.20, 1.10)
#define BORDER_COLOR vec3(1.05, 0.20, 1.00)

#define BRIGHTNESS 0.007
#define THICKNESS  0.002

//Checks if a and b are approximately equal.
bool ApproxEqual(float a, float b)
{
    return abs(a - b) <= epsilon;
}

//Distance to a line segment,
float dfLine(vec2 start, vec2 end, vec2 uv)
{
	start *= scale;
	end *= scale;

	vec2 line = end - start;
	float frac = dot(uv - start,line) / dot(line,line);
	return distance(start + line * clamp(frac, 0.0, 1.0), uv);
}

//Distance to an arc.
float dfArc(vec2 origin, float start, float sweep, float radius, vec2 uv)
{
	origin *= scale;
	radius *= scale;
	uv -= origin;

	uv *= mat2(cos(start), sin(start),-sin(start), cos(start));

    #ifdef USE_CHEAP_ARCS
        if(ApproxEqual(sweep, pi)) //180 degrees
        {
            float d = abs(length(uv) - radius) + step(uv.y, 0.0) * infinity;
            d = min(d, min(length(uv - vec2(radius, 0)), length(uv + vec2(radius, 0))));
            return d;
        }
        else if(ApproxEqual(sweep, pi/2.0)) //90 degrees
        {
            float d = abs(length(uv) - radius) + step(min(uv.x, uv.y), 0.0) * infinity;
            d = min(d, min(length(uv - vec2(0, radius)), length(uv - vec2(radius, 0))));
            return d;
        }
        else //Others
        {
            float offs = (sweep / 2.0 - pi);
            float ang = mod(atan(uv.y, uv.x) - offs, tau) + offs;
            ang = clamp(ang, min(0.0, sweep), max(0.0, sweep));

            return distance(radius * vec2(cos(ang), sin(ang)), uv);
        }
    #else
        float offs = (sweep / 2.0 - pi);
        float ang = mod(atan(uv.y, uv.x) - offs, tau) + offs;
        ang = clamp(ang, min(0.0, sweep), max(0.0, sweep));

        return distance(radius * vec2(cos(ang), sin(ang)), uv);
	#endif
}


float lower_b(vec2 uv, float dist, float offset)
{
    dist = min(dist, dfLine(vec2(0.400 + offset,1.000), vec2(0.000 + offset,0.067), uv));
    dist = min(dist, dfArc(vec2(0.033+ offset,0.333),4.992, 1.942, 0.333, uv));
    return dist;
}

float lower_l(vec2 uv, float dist, float offset)
{
    dist = min(dist, dfLine(vec2(0.400 + offset,1.000), vec2(0.000 + offset,0.067), uv));
    dist = min(dist, dfLine(vec2(0.285 + offset,1.10), vec2(0.200 + offset,0.067), uv));
    dist = min(dist, dfArc(vec2(0.337 + offset,1.05),2.742, -3.142, 0.067, uv));
    dist = min(dist, dfArc(vec2(0.287 + offset,0.05),2.742, 3.142, 0.067, uv));
    return dist;
}

float lower_p(vec2 uv, float dist, float offset)
{
    dist = min(dist, dfLine(vec2(0.000 + offset,0.400), vec2(-0.400 + offset,-0.467), uv));
    dist = min(dist, dfArc(vec2(-0.16 + offset,0.333),4.992, 1.942, 0.333, uv));

    return dist;
}

float lower_t(vec2 uv, float dist, float offset)
{
    dist = min(dist, dfLine(vec2(0.400 + offset,1.000), vec2(0.000 + offset,0.067), uv));
    dist = min(dist, dfArc(vec2(0.337 + offset,1.05),2.742, -3.142, 0.067, uv));
    dist = min(dist, dfArc(vec2(0.337 + offset,0.65),2.742, (3.142/2.0) + 0.9, 0.067, uv));
    dist = min(dist, dfArc(vec2(0.087 + offset,0.05),2.742, 3.142, 0.067, uv));


    return dist;
}

float lower_s(vec2 uv, float dist, float offset)
{
    return dist;
}

float lower_a(vec2 uv, float dist, float offset)
{

    dist = min(dist, dfArc(vec2(0.10 + offset,0.283),4.692 + pi, 2.642, 0.333, uv));
    dist = min(dist, dfLine(vec2((0.24) + offset,0.600), vec2(0.000 + offset,0.067), uv));
    dist = min(dist, dfArc(vec2(0.087 + offset,0.05),2.742, 3.142, 0.067, uv));
	dist = min(dist, dfLine(vec2((0.14) + offset,0.000), vec2(0.200 + offset,0.067), uv));
    return dist;
}

float lower_e(vec2 uv, float dist, float offset)
{
    dist = min(dist, dfLine(vec2((0.19) + offset,0.500), vec2(0.000 + offset,0.067), uv));
    dist = min(dist, dfArc(vec2(0.40 + offset,0.5),2.942, -(pi-1.3), 0.187, uv));
    dist = min(dist, dfLine(vec2((0.27) + offset,0.300), vec2(0.3 + offset,0.367), uv));
    dist = min(dist, dfArc(vec2(0.24 + offset,0.09),3.6, (1.5), 0.087, uv));
    dist = min(dist, dfLine(vec2(0.370 + offset,0.067), vec2(0.31 + offset,0.000), uv));
    return dist;
}

float lower_i(vec2 uv, float dist, float offset)
{
    dist = min(dist, dfLine(vec2((0.24) + offset,0.600), vec2(0.000 + offset,0.067), uv));
    dist = min(dist, dfArc(vec2(0.087 + offset,0.05),2.742, 3.142, 0.067, uv));
    dist = min(dist, dfLine(vec2((0.34) + offset,0.700), vec2(0.34 + offset,0.797), uv));
	dist = min(dist, dfLine(vec2((0.14) + offset,0.000), vec2(0.200 + offset,0.067), uv));

    return dist;
}

float lower_u(vec2 uv, float dist, float offset)
{
    dist = min(dist, dfLine(vec2((0.24) + offset,0.600), vec2(0.000 + offset,0.067), uv));
    dist = min(dist, dfArc(vec2(0.12 + offset,0.05),2.842, 2.942, 0.097, uv));
	dist = min(dist, dfLine(vec2((0.24) + offset +0.23,0.600), vec2(0.230 + offset,0.067), uv));
	dist = min(dist, dfArc(vec2(0.287 + offset,0.05),2.742, 3.142, 0.067, uv));
    dist = min(dist, dfLine(vec2((0.34) + offset,0.000), vec2(0.400 + offset,0.067), uv));

    return dist;
}



float lower_c(vec2 uv, float dist, float offset)
{
    dist = min(dist, dfLine(vec2((0.19) + offset,0.500), vec2(0.000 + offset,0.067), uv));
    dist = min(dist, dfArc(vec2(0.40 + offset,0.5),2.942, -(pi-1.3), 0.187, uv));
    //dist = min(dist, dfLine(vec2((0.27) + offset,0.300), vec2(0.3 + offset,0.367), uv));
    dist = min(dist, dfArc(vec2(0.24 + offset,0.09),3.6, (1.5), 0.087, uv));
    dist = min(dist, dfLine(vec2(0.370 + offset,0.067), vec2(0.31 + offset,0.000), uv));
    return dist;
}


float upper_s(vec2 uv, float dist, float offset)
{
    dist = min(dist, dfLine(vec2(0.267 + offset,1.200), vec2(0.533 + offset,1.200), uv));
	dist = min(dist, dfLine(vec2(0.267 + offset,0.667), vec2(0.533 + offset,0.667), uv));
	dist = min(dist, dfLine(vec2(0.533 + offset,0.000), vec2(0.067 + offset,0.000), uv));
	dist = min(dist, dfLine(vec2(0.400 + offset,0.133), vec2(0.067 + offset,0.133), uv));
    dist = min(dist, dfArc(vec2(0.267 + offset,0.933),1.571, 3.142, 0.267, uv));
	dist = min(dist, dfArc(vec2(0.067 + offset,0.067),1.571, 3.142, 0.067, uv));
	dist = min(dist, dfArc(vec2(0.533 + offset,0.333),4.712, 3.142, 0.333, uv));

	return dist;
}

float dfLogo(vec2 uv)
{
	float dist = infinity;

 	dist = lower_b(uv, dist, 0.0);
  dist = lower_l(uv, dist, 0.4);
  dist = lower_u(uv, dist, 0.8);
  dist = lower_e(uv, dist, 1.2);

  dist = lower_p(uv, dist, 2.2);
  dist = lower_l(uv, dist, 2.4);
  dist = lower_a(uv, dist, 3.0);
  dist = lower_t(uv, dist, 3.2);
  dist = lower_e(uv, dist, 3.6);

  dist = upper_s(uv, dist, 4.5);
  dist = lower_p(uv, dist, 5.7);
  dist = lower_e(uv, dist, 5.9);
  dist = lower_c(uv, dist, 6.3);
  dist = lower_i(uv, dist, 6.7);
  dist = lower_a(uv, dist, 7.2);
  dist = lower_l(uv, dist, 7.4);

	return dist;
}

float dfBorder(vec2 uv)
{
    float dist = infinity;

	dist = min(dist, dfLine(vec2(0.067,1.533), vec2(8.733,1.533), uv));
	dist = min(dist, dfLine(vec2(9.133,1.133), vec2(9.133,0.067), uv));
	dist = min(dist, dfLine(vec2(8.733,-0.333), vec2(4.467,-0.333), uv));
	dist = min(dist, dfLine(vec2(-0.333,0.067), vec2(-0.333,1.133), uv));
	dist = min(dist, dfLine(vec2(0.067,1.400), vec2(4.333,1.400), uv));
	dist = min(dist, dfLine(vec2(9.000,1.133), vec2(9.000,0.067), uv));
	dist = min(dist, dfLine(vec2(8.733,-0.200), vec2(0.067,-0.200), uv));
	dist = min(dist, dfLine(vec2(-0.200,0.067), vec2(-0.200,1.133), uv));
	dist = min(dist, dfLine(vec2(4.333,-0.333), vec2(0.067,-0.333), uv));
	dist = min(dist, dfLine(vec2(4.467,1.400), vec2(8.733,1.400), uv));
	dist = min(dist, dfArc(vec2(8.733,1.133),0.000, 1.571, 0.400, uv));
	dist = min(dist, dfArc(vec2(8.733,0.067),4.712, 1.571, 0.400, uv));
	dist = min(dist, dfArc(vec2(0.067,0.067),3.142, 1.571, 0.400, uv));
	dist = min(dist, dfArc(vec2(0.067,1.133),1.571, 1.571, 0.400, uv));
	dist = min(dist, dfArc(vec2(8.733,1.133),0.000, 1.571, 0.267, uv));
	dist = min(dist, dfArc(vec2(8.733,0.067),4.712, 1.571, 0.267, uv));
	dist = min(dist, dfArc(vec2(0.067,0.067),3.142, 1.571, 0.267, uv));
	dist = min(dist, dfArc(vec2(0.067,1.133),1.571, 1.571, 0.267, uv));

    return dist;
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 aspect = iResolution.xy / iResolution.y;
	vec2 uv = fragCoord.xy / iResolution.y - aspect/2.0;

    vec2 offs = vec2(9.0, 1.0) * scale/2.0;

    float dist = 0.0;
    float shade = 0.0;
    vec3 color = vec3(0);

    //Flicker fade in effect.
    float tf_text = max(epsilon, iTime - 0.6);
    float bright_text = BRIGHTNESS * min(1.0, 1.0 - sin(tf_text * pi * 50.0) / (tf_text * pi * 1.3));

    float tf_bord = max(epsilon, iTime - 0.5);
    float bright_bord = BRIGHTNESS * min(1.0, 1.0 - sin(tf_bord * pi * 50.0) / (tf_bord * pi * 1.3));

    //"Shadertoy"
	dist = dfLogo(uv + offs);

	shade = bright_text / max(epsilon, dist - THICKNESS);

	color += TEXT_COLOR * shade;

    //Border
    dist = dfBorder(uv + offs);

	shade = bright_bord / max(epsilon, dist - THICKNESS);

	color += BORDER_COLOR * shade;
	fragColor = vec4(color , 1.0);
}
