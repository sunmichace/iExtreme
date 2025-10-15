
function draw_sunbrust(data){
console.log("sunbrust");
var statistics_chart = echarts.init(document.getElementById('sunburst-frame'));
var option;
// var data = ;

option = {
    title: {
      text: 'Virus Factor Analysis',
      left: "center",
      subtext: '',
      textStyle: {
        fontSize: 36,
        align: 'center'
      },
      subtextStyle: {
        align: 'center'
      },
      sublink: ''
    },
    tooltip: {
        trigger: 'item',
        formatter: function(params) {  
            //console.log(params.value);
            var temp = `<span style='font-weight:bold;'>${params.name}</span> Count: ${params.value}`;
            return temp;  
        },
    },
    series: {
      type: 'sunburst',
      data: data,
      radius: [0, '95%'],
      sort: undefined,
      colorBy: "data",
      emphasis: {
        focus: 'ancestor'
      },
      levels: [
        {},
        {
          r0: '15%',
          r: '35%',
          itemStyle: {
            borderWidth: 2
          },
          label: {
            rotate: 'tangential',
            fontSize: 18
          }
        },
        
        {
          r0: '35%',
          r: '70%',
          label: {
            position: 'outside',
            padding: -50,
            silent: false 
          },
          itemStyle: {
            borderWidth: 3
          }
        }
      ]
    }
};
statistics_chart.setOption(option);
}
