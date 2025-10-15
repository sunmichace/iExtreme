
function draw_TSpH(temp_x,temp_y,salinity_x,salinity_y,pH_x,pH_y){
  var statistics_chart = echarts.init(document.getElementById('TSpH-frame'));
  var option;

  const colors = ['#E4080A', '#1E90FF', '#7DDA58'];
option = {
  color: colors,
  tooltip: {
    trigger: 'none',
    axisPointer: {
      type: 'cross'
    }
  },
  legend: {},
  grid: {
    top: 70,
    bottom: 120  // 增加底部空间
  },
  xAxis: [
    {
      type: 'category',
      position: 'bottom',
      offset: 0,
      axisLabel: {
        fontSize: 20 
      },
      axisTick: {
        alignWithLabel: true
      },
      axisLine: {
        onZero: false,
        lineStyle: {
          color: colors[0]
        }
      },
      axisPointer: {
        label: {
          formatter: function (params) {
            return (
              'Temperature  ' +
              params.value +
              (params.seriesData.length ? ' Confidence: ' + params.seriesData[0].data : '')
            );
          }
        }
      },
      // prettier-ignore
      data: temp_x
    },
    {
      type: 'category',
      position: 'bottom',
      offset: 40,
      axisLabel: {
        fontSize: 20 
      },
      axisTick: {
        alignWithLabel: true
      },
      axisLine: {
        onZero: false,
        lineStyle: {
          color: colors[1]
        }
      },
      axisPointer: {
        label: {
          formatter: function (params) {
            return (
              'Salinity  ' +
              params.value +
              (params.seriesData.length ? ' Confidence: ' + params.seriesData[0].data : '')
            );
          }
        }
      },
      // prettier-ignore
      data: salinity_x
    },
    {
      type: 'category',
      position: 'bottom',
      offset: 80,
      axisLabel: {
        fontSize: 20 
      },
      axisTick: {
        alignWithLabel: true
      },
      axisLine: {
        onZero: false,
        lineStyle: {
          color: colors[2]
        }
      },
      axisPointer: {
        label: {
          formatter: function (params) {
            return (
              'pH  ' +
              params.value +
              (params.seriesData.length ? ' Confidence: ' + params.seriesData[0].data : '')
            );
          }
        }
      },
      // prettier-ignore
      data: pH_x
    }
  ],
  yAxis: [
    {
      name: 'Confidence(%)',
      nameLocation: 'end',
      nameTextStyle: {
        fontSize: 20
      },
      axisLabel: {
        fontSize: 20 
      },
    }
  ],
  series: [
    {
      name: 'Temperature',
      type: 'line',
      xAxisIndex: 0,
      smooth: true,
      lineStyle:{
        width: 5
      },
      emphasis: {
        focus: 'series'
      },
      data: temp_y
    },
    {
      name: 'Salinity',
      type: 'line',
      xAxisIndex: 1,
      smooth: true,
      lineStyle:{
        width: 5
      },
      emphasis: {
        focus: 'series'
      },
      data: salinity_y
    },
    {
      name: 'pH',
      type: 'line',
      xAxisIndex: 2,
      smooth: true,
      lineStyle:{
        width: 5
      },
      emphasis: {
        focus: 'series'
      },
      data: pH_y
    }
  ]
};
  statistics_chart.setOption(option);
}

// T_x = [10,63.05,67.12,71.2,75.27,79.34,83.42,87.49,91.57,95.64]
// T_y = [0.25,0.38,0.75,2.77,12.2,82.01,17.99,5.66,1.13,0.38]
// S_x = [0.,0.,1.36,3.5,5.65,7.8,9.95,12.09,14.24,16.39]
// S_y = [0.54,3.15,10.11,89.89,7.39,3.15,1.2,0.87,0.43,0.22]
// pH_x = [3.9,4.19,4.49,4.78,5.08,5.37,5.67,5.96,6.26,6.55]
// pH_y = [1.04, 2.43, 6.6,12.5,50.,50.,10.42,5.56,3.12,0.69]
// draw_TSpH(T_x, T_y, S_x, S_y, pH_x, pH_y);
// var test_data = [0.12,0.22,0.95]
// draw_leida(test_data)