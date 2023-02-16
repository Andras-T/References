package hu.bme.aut.steptracker.fragments

import android.graphics.Color
import android.os.Bundle
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import androidx.fragment.app.Fragment
import com.github.aachartmodel.aainfographics.aachartcreator.AAChartAnimationType
import com.github.aachartmodel.aainfographics.aachartcreator.AAChartModel
import com.github.aachartmodel.aainfographics.aachartcreator.AAChartType
import com.github.aachartmodel.aainfographics.aachartcreator.AASeriesElement
import hu.bme.aut.steptracker.databinding.FragmentDataBinding
import hu.bme.aut.steptracker.model.Score
import hu.bme.aut.steptracker.sqlite.PersistentDataHelper
import java.util.*


class DataFragment : Fragment() {

    private lateinit var binding : FragmentDataBinding
    private lateinit var values : Array<Any>
    private lateinit var dataHelper: PersistentDataHelper
    private val weekDays: Array<String> = arrayOf("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

    override fun onCreateView(inflater: LayoutInflater, container: ViewGroup?, savedInstanceState: Bundle?): View {
        super.onCreateView(inflater, container, savedInstanceState)
        binding = FragmentDataBinding.inflate(layoutInflater, container, false)
        dataHelper = activity?.let { PersistentDataHelper(it) }!!
        return binding.root
    }

    override fun onViewCreated(view: View, savedInstanceState: Bundle?) {
        super.onViewCreated(view, savedInstanceState)
        values = arrayOf(0,0,0,0,0,0,0)
        readData()
        drawChart()
    }

    private fun readData() {
        dataHelper.open()
        val scores = dataHelper.searchWeeklyScore(Calendar.getInstance().get(Calendar.WEEK_OF_YEAR).toString())
        for (score: Score in scores) {
            //TimeFormat: Sat Nov 19 2022 47
            val day = score.date.dropLast(15)
            val index = weekDays.indexOf(day)
            values[index] = score.steps
        }
        dataHelper.close()
    }

    private fun drawChart() {
        val aaChartView = binding.aaChartView
        val aaChartModel : AAChartModel = AAChartModel()
            .chartType(AAChartType.Column)
            .animationType(AAChartAnimationType.Bounce)
            .title("Steps taken each day")
            .subtitle("Current Week")
            .tooltipValueSuffix(" steps")
            .backgroundColor(Color.CYAN)
            .categories(weekDays)
            .series(arrayOf(
                AASeriesElement()
                    .name("Steps taken")
                    .data(values)))
        aaChartView.aa_drawChartWithChartModel(aaChartModel)
        aaChartView.invalidate()
    }
}