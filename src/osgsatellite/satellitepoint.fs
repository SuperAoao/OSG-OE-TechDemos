// version OpenGL 2.1 & GLSL 1.2, 
varying float window_pos_x;
varying float window_pos_y;

uniform float n_pe;
varying float ptSize;
varying float highlightType;
uniform float n_e;  // 每1个灰度值对应的电子数
uniform sampler2D tex;
void main()
{
    // gl_FragCoord holds the window relative coordinates x, y, z, and 1/w values for the fragment
    gl_FragColor.a = 1.0;
    gl_FragColor.r = 0.0;
    if (window_pos_x - gl_FragCoord.x <1e-3)
    {
        gl_FragColor.r = 1.0;
    }
  
    gl_FragColor.g = 1.0;
    float temp = texture(tex, gl_PointCoord).x;
    //temp = n_pe * temp / n_e / 256.0; 
    gl_FragColor = vec4(temp, temp, temp, 1.0);
    // if (temp < 0.02)
    // {
    //     gl_FragColor.a = 0.0;
    // }
}
