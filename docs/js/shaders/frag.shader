//==============================================================================
// Base P5 Fragment Shader
//==============================================================================
#ifdef GL_ES
precision mediump float; // set precision level if available
#endif
//==============================================================================
// Global Constants
const float twoPi = 6.283185307179586;
//==============================================================================
// input variables
uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;
//==============================================================================
// varying: these have all come from the vertex shader
varying vec3 var_vertPos;
varying vec4 var_vertCol;
varying vec3 var_vertNormal;
varying vec2 var_vertTexCoord;
varying vec4 v_color;
//==============================================================================

const float maxz = 1.0/50.0;
void main( void )
{
    float color = var_vertPos.z*maxz;
    gl_FragColor = vec4(
    clamp(0.55 * color,0.,1.) + clamp(0.7*-color,0.,1.) + 0.17,
    clamp(0.35 * color,0.,1.) + clamp(0.2*-color,0.,1.) + 0.17,
    clamp(0.25 * color,0.,1.) + clamp(0.7*-color,0.,1.) + 0.11,
    1.0);
}
