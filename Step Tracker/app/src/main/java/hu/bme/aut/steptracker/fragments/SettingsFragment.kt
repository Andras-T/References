package hu.bme.aut.steptracker.fragments

import android.content.Context
import android.os.Bundle
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import androidx.appcompat.app.AppCompatDelegate
import androidx.appcompat.app.AppCompatDelegate.MODE_NIGHT_NO
import androidx.appcompat.app.AppCompatDelegate.MODE_NIGHT_YES
import androidx.fragment.app.Fragment
import com.google.android.material.snackbar.Snackbar
import hu.bme.aut.steptracker.databinding.FragmentSettingsBinding

class SettingsFragment : Fragment() {

    private lateinit var binding : FragmentSettingsBinding

    override fun onCreateView(inflater: LayoutInflater, container: ViewGroup?, savedInstanceState: Bundle?): View {
        super.onCreateView(inflater, container, savedInstanceState)
        binding = FragmentSettingsBinding.inflate(inflater, container, false)
        applyListener()
        return binding.root
    }

    private fun applyListener() {
        val sharedPreference =  activity?.getSharedPreferences("Theme", Context.MODE_PRIVATE)
        val storedTheme: String? = sharedPreference?.getString("Theme","defaultName")
        if (storedTheme!=null && storedTheme =="Night"){
            binding.changeTheme.isChecked = true
        }
        var theme: String
        binding.changeTheme.setOnCheckedChangeListener { buttonView, isChecked ->
            theme = if (isChecked){
                AppCompatDelegate.setDefaultNightMode(MODE_NIGHT_YES)
                "Night"
            } else {
                AppCompatDelegate.setDefaultNightMode(MODE_NIGHT_NO)
                "Normal"
            }

            val editor = sharedPreference?.edit()
            editor?.putString("Theme",theme)
            editor?.apply()
        }

        binding.submitBtn.setOnClickListener {
            val height = binding.heightInput.text.toString()
            val weight = binding.weightInput.text.toString()
            if (height.isNotEmpty() && weight.isNotEmpty())
            {
                if (height.toInt() < 50 || weight.toInt() < 50)
                {
                    Snackbar.make(binding.root, "Minimum values are 50", Snackbar.LENGTH_SHORT).show()
                }
                else {
                    submitData(height, weight)
                }
            }
            else {
                Snackbar.make(binding.root, "Missing data", Snackbar.LENGTH_SHORT).show()
            }

        }

    }
    private fun submitData(height: String, weight: String) {
        val sharedPreferenceHeight =  activity?.getSharedPreferences("Height", Context.MODE_PRIVATE)
        val editorHeight = sharedPreferenceHeight?.edit()
        editorHeight?.putString("Height", height)
        editorHeight?.apply()
        val sharedPreferenceWeight =  activity?.getSharedPreferences("Weight", Context.MODE_PRIVATE)
        val editorWeight = sharedPreferenceWeight?.edit()
        editorWeight?.putString("Weight", weight)
        editorWeight?.apply()
        Snackbar.make(binding.root, "Submited", Snackbar.LENGTH_SHORT).show()
    }

}