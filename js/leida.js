
function draw_leida(data){
  var statistics_chart = echarts.init(document.getElementById('leida-frame'));
  var option;

  option = {
  title: {
    text: 'iExtreme prediction result',
    left: "center",
    textStyle: {
      fontSize: 32,
      align: 'center'
    },
  },
  
  radar: {
    shape: '',
    indicator: [
      { name: 'Thermophilic (>70℃)', max: 1 },
      { name: 'Halophilic (>15%)', max: 1 },
      { name: 'Acidophilic (pH<5)', max: 1 }
    ],
    nameGap:30,
    axisName: {
      //formatter: '【{value}】',
      color: '#fff',
      backgroundColor: '#CECECE',
      borderRadius: 3,
      padding: [3, 5],
      fontSize: 24,
    }
  },
  series: [
    {
      name: 'Score',
      type: 'radar',
      data: [
        {
          value: data,
          name: 'user',
          lineStyle: {
            type: 'dashed'
          },
          label: {
            show: true,
            fontSize: 18,
            position: 'left',
            color: 'black',
            fontStyle: 'bold',
            formatter: function (params) {
              return params.value;
            }
          },
          areaStyle: {
            color: 'rgba(255, 145, 124, 0.9)'
          }
        }
      ]
    }
  ]
  };
  statistics_chart.setOption(option);
}

// var test_data = [0.12,0.22,0.95]
// draw_leida(test_data)